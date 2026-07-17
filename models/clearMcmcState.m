function clearMcmcState(storageDir, modelName)

path = fullfile(storageDir, sprintf('%s_mcmcState.mat', modelName));
if isfile(path)
  delete(path);
end
legacy = fullfile(storageDir, sprintf('%s_thinState.mat', modelName));
if isfile(legacy)
  delete(legacy);
end
end
