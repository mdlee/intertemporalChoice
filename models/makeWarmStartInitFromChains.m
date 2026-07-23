function initGenerator = makeWarmStartInitFromChains(chains, monitorParams, nParticipants, varargin)
%MAKEWARMSTARTINITFROMCHAINS  Init generator near posterior means of a prior fit.
%
%   INITGENERATOR = makeWarmStartInitFromChains(CHAINS, MONITORPARAMS, NPARTICIPANTS)
%   returns a function handle suitable for trinity/callbayes 'init'. Each call
%   (one per chain) returns participant-level monitored parameters jittered
%   around the posterior means in CHAINS. Values are clamped to the
%   truncated-normal supports used in the hierarchical entrop JAGS models.
%
%   Note: group-level mu_*/tau_* are not initialized here. Trinity's JAGS
%   writer turns underscores into dots in .init files, so those names would
%   not match models that use allowunderscores=1.
%
%   Name-value:
%     jitterFrac (0.05) — relative Gaussian jitter scale (per chain)
%     seed      ([])    — optional RNG seed for reproducibility of first draws

p = inputParser;
addParameter(p, 'jitterFrac', 0.05, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
parse(p, varargin{:});
jitterFrac = p.Results.jitterFrac;

monitorParams = cellstr(monitorParams(:));
nParticipants = double(nParticipants);
means = struct();
for k = 1:numel(monitorParams)
  pname = monitorParams{k};
  mu = get_matrix_from_coda(chains, pname, @mean);
  mu = mu(:);
  if numel(mu) < nParticipants
    error('makeWarmStartInitFromChains:badChains', ...
      'Expected %d elements for %s, found %d.', nParticipants, pname, numel(mu));
  end
  means.(pname) = mu(1:nParticipants);
end

if ~isempty(p.Results.seed)
  rng(p.Results.seed);
end

initGenerator = @() drawWarmStart(means, monitorParams, jitterFrac);
end

function s = drawWarmStart(means, monitorParams, jitterFrac)
s = struct();
for k = 1:numel(monitorParams)
  pname = monitorParams{k};
  base = means.(pname);
  [lo, hi] = paramBounds(pname);
  jittered = base + jitterFrac .* max(abs(base), 0.05) .* randn(size(base));
  s.(pname) = clampVec(jittered, lo, hi);
end
end

function [lo, hi] = paramBounds(pname)
switch pname
  case {'gamma', 'tau', 'kappa', 'vartheta'}
    lo = 1e-4; hi = 100;
  case 'eta'
    lo = 0; hi = 500;
  case 'w'
    lo = 0; hi = 20;
  case {'beta0', 'betaRA', 'betaRR', 'betaTA', 'betaTR'}
    lo = -2; hi = 2;
  otherwise
    lo = -inf; hi = inf;
end
end

function x = clampVec(x, lo, hi)
x = min(max(x, lo), hi);
% Keep strictly inside open truncations used by JAGS T(,)
epsIn = 1e-6;
if isfinite(lo)
  x = max(x, lo + epsIn);
end
if isfinite(hi)
  x = min(x, hi - epsIn);
end
x = x(:);
end
