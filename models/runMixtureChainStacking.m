%% runMixtureChainStacking — Yao/Vehtari/Gelman (2022) stacking on mixture chains
%
% The sequential mixture run only stacks z for R-hat. Stacking needs full θ
% draws, so this script:
%   1. Resumes last-state inits from the sequential checkpoint
%   2. Collects a short full-parameter sample per chain (does not touch the
%      sequential stack / stop file)
%   3. Optional: cluster chains by categorical z occupancy
%   4. PSIS-LOO pointwise predictive densities per chain/cluster
%   5. Stacking weights maximizing LOO log score (Yao et al. 2022)
%
%   runMixtureChainStacking

clear; close all;

modelsDir = fileparts(mfilename('fullpath'));
cd(modelsDir);
addpath(modelsDir);
addpath(fullfile(modelsDir, '..', 'general'));

logDir = fullfile(modelsDir, 'logs');
storageDir = fullfile(modelsDir, 'storage');
figuresDir = fullfile(modelsDir, 'figures');
if ~isfolder(figuresDir), mkdir(figuresDir); end

dataName = 'intertemporalChoice';
engine = 'jags';
modelName = 'latentMixtureHierarchicalPrecision_entrop';
checkpointFile = fullfile(logDir, 'runLatentMixtureSequential.checkpoint.mat');
outFile = fullfile(storageDir, 'latentMixtureChainStacking.mat');
nModels = 11;
modelNames = {'EX','HC','HY','PD','DD','TR','UT','ITCH','Guess','LL','SS'};

nSamples = 1500;   % draws per chain for stacking (full θ)
nBurnin = 200;
nThin = 1;
doCluster = true;
rhatCluster = 1.1; % merge chains if pairwise z R-hat of modes agree

fprintf('=== Mixture chain stacking (Yao, Vehtari & Gelman 2022) ===\n');
if ~isfile(checkpointFile)
  error('Missing checkpoint: %s', checkpointFile);
end
ck = load(checkpointFile);
fprintf('Checkpoint: stacked z n=%d batch=%d best z R-hat=%.3f\n', ...
  ck.nCollected, ck.batch, ck.bestRhat);
if ~isfield(ck, 'lastInitCell') || isempty(ck.lastInitCell)
  error('Checkpoint has no lastInitCell — cannot continue chains.');
end
nChains = numel(ck.lastInitCell);
initFn = makeSerialInitFnLocal(ck.lastInitCell);

[data, d] = prepareIntertemporalChoiceData(dataName);
data.nModels = nModels;

% Monitor the same participant + group fields as the sequential runner
s0 = ck.lastInitCell{1};
monitorParams = fieldnames(s0);
fprintf('Collecting %d draws × %d chains (burn-in %d) with full monitors...\n', ...
  nSamples, nChains, nBurnin);

tic;
modelPath = resolveJagsModelFile(modelsDir, sprintf('%s_%s.txt', modelName, engine));
[~, chains, ~, ~] = callbayes(engine, ...
  'model', modelPath, ...
  'data', data, ...
  'outputname', 'samples', ...
  'init', initFn, ...
  'datafilename', [modelName '_stack'], ...
  'initfilename', [modelName '_stack'], ...
  'scriptfilename', [modelName '_stack'], ...
  'logfilename', fullfile('tmp', [modelName '_stack']), ...
  'nchains', nChains, ...
  'nburnin', nBurnin, ...
  'nsamples', nSamples, ...
  'monitorparams', monitorParams, ...
  'thin', nThin, ...
  'workingdir', fullfile('tmp', [modelName '_stack']), ...
  'verbosity', 0, ...
  'saveoutput', true, ...
  'allowunderscores', 1, ...
  'parallel', true);
fprintf('JAGS took %.1f s\n', toc);

% --- optional clustering by z occupancy ---
clusterId = (1:nChains)';
if doCluster
  occ = zeros(nChains, nModels);
  for c = 1:nChains
    for pp = 1:d.nParticipants
      zc = round(chains.(sprintf('z_%d', pp))(:, c));
      for m = 1:nModels
        occ(c, m) = occ(c, m) + mean(zc == m);
      end
    end
    occ(c, :) = occ(c, :) / d.nParticipants;
  end
  % Merge chains with total-variation distance < 0.15 in mean occupancy
  unused = true(nChains, 1);
  cid = 0;
  clusterId = zeros(nChains, 1);
  for c = 1:nChains
    if ~unused(c), continue; end
    cid = cid + 1;
    clusterId(c) = cid;
    unused(c) = false;
    for c2 = c + 1:nChains
      if ~unused(c2), continue; end
      tv = 0.5 * sum(abs(occ(c, :) - occ(c2, :)));
      if tv < 0.15
        clusterId(c2) = cid;
        unused(c2) = false;
      end
    end
  end
  fprintf('Clustering: %d chains → %d clusters (TV<0.15 on mean z occupancy)\n', ...
    nChains, max(clusterId));
  for k = 1:max(clusterId)
    fprintf('  cluster %d: chains [%s]\n', k, num2str(find(clusterId == k)'));
  end
else
  fprintf('Clustering skipped — one component per chain\n');
end
K = max(clusterId);

% --- PSIS-LOO LPD per observation × cluster ---
nObs = d.nParticipants * d.nTrials;
looLpd = nan(nObs, K);
khatMax = nan(1, K);
fprintf('Computing PSIS-LOO LPD (%d obs × %d clusters)...\n', nObs, K);
tic;
for k = 1:K
  members = find(clusterId == k);
  % Pool draws across member chains
  chK = subsetChainCols(chains, members);
  % Flatten to single long chain for LOO (nDraws pooled)
  chFlat = flattenChains(chK);
  nSk = size(chFlat.z_1, 1);
  fprintf('  cluster %d: %d chains, %d pooled draws\n', k, numel(members), nSk);

  ll = latentMixturePointwiseLogLik(chFlat, data);  % nObs × nSk
  for i = 1:nObs
    [logW, kh] = psisSmoothLogWeights(-ll(i, :)');  % ratios 1/p
    khatMax(k) = max(khatMax(k), kh);
    % pk(yi|y-i) ≈ sum p * r / sum r  with r = exp(logW)
    logP = ll(i, :)';
    logNum = logsumexpVec(logP + logW);
    logDen = logsumexpVec(logW);
    looLpd(i, k) = logNum - logDen;
  end
  fprintf('    max khat=%.3f  mean loo lpd=%.4f\n', khatMax(k), mean(looLpd(:, k)));
end
fprintf('PSIS-LOO took %.1f s\n', toc);

% --- stacking weights ---
[wStack, elpdStack, info] = stackingOptimizeWeights(looLpd, 'dirichletAlpha', 1.01);
wUnif = ones(1, K) / K;
elpdUnif = info.uniformElpd;

fprintf('\n=== Stacking weights ===\n');
for k = 1:K
  fprintf('  cluster %d (chains [%s]): w=%.4f\n', k, num2str(find(clusterId == k)'), wStack(k));
end
fprintf('ELPD stacking = %.2f\n', elpdStack);
fprintf('ELPD uniform  = %.2f   (Δ = %.2f)\n', elpdUnif, elpdStack - elpdUnif);
fprintf('PSIS max khat by cluster: %s\n', mat2str(khatMax, 3));

% --- stacked vs equal-weight z posterior (from this sample) ---
zPostStack = zeros(d.nParticipants, nModels);
zPostUnif = zeros(d.nParticipants, nModels);
for pp = 1:d.nParticipants
  for k = 1:K
    members = find(clusterId == k);
    zc = [];
    for c = members(:)'
      zc = [zc; round(chains.(sprintf('z_%d', pp))(:, c))]; %#ok<AGROW>
    end
    zc = max(1, min(nModels, zc));
    pk = accumarray(zc, 1, [nModels, 1]) / numel(zc);
    zPostStack(pp, :) = zPostStack(pp, :) + wStack(k) * pk(:)';
    zPostUnif(pp, :) = zPostUnif(pp, :) + wUnif(k) * pk(:)';
  end
end

fprintf('\nMAP model under stacking vs uniform (participants with disagreement):\n');
for pp = 1:d.nParticipants
  [~, ms] = max(zPostStack(pp, :));
  [~, mu] = max(zPostUnif(pp, :));
  if ms ~= mu || pp == 11 || pp == 19
    fprintf('  p%02d  stack=%s (%.2f)  unif=%s (%.2f)\n', pp, ...
      modelNames{ms}, zPostStack(pp, ms), modelNames{mu}, zPostUnif(pp, mu));
  end
end

% Highlight p11 / p19
fprintf('\np11 model probs (stack / unif):\n');
for m = 1:nModels
  if zPostStack(11, m) > 0.01 || zPostUnif(11, m) > 0.01
    fprintf('  %s  %.3f / %.3f\n', modelNames{m}, zPostStack(11, m), zPostUnif(11, m));
  end
end
fprintf('p19 model probs (stack / unif):\n');
for m = 1:nModels
  if zPostStack(19, m) > 0.01 || zPostUnif(19, m) > 0.01
    fprintf('  %s  %.3f / %.3f\n', modelNames{m}, zPostStack(19, m), zPostUnif(19, m));
  end
end

save(outFile, 'wStack', 'wUnif', 'elpdStack', 'elpdUnif', 'looLpd', ...
  'clusterId', 'khatMax', 'zPostStack', 'zPostUnif', 'chains', ...
  'nSamples', 'nBurnin', 'modelNames', 'info', 'ck', '-v7.3');
fprintf('\nSaved %s\n', outFile);
fprintf('=== done ===\n');

function fh = makeSerialInitFnLocal(initCell)
idx = 0;
fh = @nextInit;
  function s = nextInit
    idx = idx + 1;
    if idx > numel(initCell)
      error('makeSerialInitFnLocal:exhausted', 'Init list exhausted at %d.', idx);
    end
    s = initCell{idx};
  end
end

function out = subsetChainCols(chains, cols)
fn = fieldnames(chains);
out = struct();
for k = 1:numel(fn)
  x = chains.(fn{k});
  if isnumeric(x) && size(x, 2) >= max(cols)
    out.(fn{k}) = x(:, cols);
  else
    out.(fn{k}) = x;
  end
end
end

function out = flattenChains(chains)
% nS × nC → (nS*nC) × 1
fn = fieldnames(chains);
out = struct();
for k = 1:numel(fn)
  x = chains.(fn{k});
  if isnumeric(x) && ndims(x) == 2
    out.(fn{k}) = x(:);
  else
    out.(fn{k}) = x;
  end
end
end

function y = logsumexpVec(x)
x = x(:);
m = max(x);
y = m + log(sum(exp(x - m)));
end
