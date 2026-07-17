function [nThin, nBurnin] = initialMcmcFromState(storageDir, modelName, defaultThin, defaultBurnin)
%INITIALMCMCFROMSTATE  Resume after failure: thin ×2; burn-in unchanged.

nThin = defaultThin;
nBurnin = defaultBurnin;
path = mcmcStatePath(storageDir, modelName);
if ~isfile(path)
  legacy = fullfile(storageDir, sprintf('%s_thinState.mat', modelName));
  if isfile(legacy)
    path = legacy;
  else
    return;
  end
end
S = load(path);
if isfield(S, 'lastUnsuccessfulThin') && isfinite(S.lastUnsuccessfulThin) && S.lastUnsuccessfulThin >= 1
  nThin = nextThin(S.lastUnsuccessfulThin);
end
if nThin ~= defaultThin
  fprintf('%s: resuming thin=%d, burn-in=%g\n', modelName, nThin, nBurnin);
end
end

function path = mcmcStatePath(storageDir, modelName)
path = fullfile(storageDir, sprintf('%s_mcmcState.mat', modelName));
end
