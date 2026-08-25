function [w, elpd, info] = stackingOptimizeWeights(looLpd, varargin)
%STACKINGOPTIMIZEWEIGHTS  Yao/Vehtari/Gelman stacking weights from LOO LPD.
%
%   LOOLPD is nObs-by-K matrix of leave-one-out log predictive densities.
%   Maximizes sum_i log(sum_k w_k exp(looLpd(i,k))) over the simplex,
%   optionally with a Dirichlet prior.

p = inputParser;
addParameter(p, 'dirichletAlpha', 1.01, @isnumeric);
addParameter(p, 'w0', [], @(x) isempty(x) || isnumeric(x));
parse(p, varargin{:});
alpha0 = p.Results.dirichletAlpha;

[nObs, K] = size(looLpd);
if K < 2
  w = 1;
  elpd = sum(looLpd);
  info = struct('exitflag', 1, 'K', K);
  return;
end

% Unconstrained α in R^{K-1}; w = softmax([α; 0])
if isempty(p.Results.w0)
  a0 = zeros(K - 1, 1);
else
  w0 = p.Results.w0(:)' ;
  w0 = w0 / sum(w0);
  a0 = log(max(w0(1:K-1), 1e-12)) - log(max(w0(K), 1e-12));
  a0 = a0(:);
end

opts = optimoptions('fminunc', ...
  'Display', 'off', ...
  'Algorithm', 'quasi-newton', ...
  'MaxFunctionEvaluations', 5000);

obj = @(a) stackingNegElpd(a, looLpd, alpha0);
[aHat, fval, exitflag] = fminunc(obj, a0, opts);
w = softmaxStack([aHat(:); 0]);
elpd = -fval;
% Remove Dirichlet contribution from reported elpd for interpretability
elpd = sum(logsumexpRows(looLpd + log(w(:)')));

info = struct('exitflag', exitflag, 'K', K, 'nObs', nObs, ...
  'dirichletAlpha', alpha0, 'uniformElpd', sum(logsumexpRows(looLpd + log(1/K))));
end

function neg = stackingNegElpd(a, looLpd, alpha0)
K = size(looLpd, 2);
w = softmaxStack([a(:); 0]);
lpd = logsumexpRows(looLpd + log(w(:)'));
prior = (alpha0 - 1) * sum(log(max(w, realmin)));
neg = -(sum(lpd) + prior);
end

function w = softmaxStack(x)
x = x - max(x);
ex = exp(x);
w = ex / sum(ex);
end

function y = logsumexpRows(X)
m = max(X, [], 2);
y = m + log(sum(exp(X - m), 2));
end
