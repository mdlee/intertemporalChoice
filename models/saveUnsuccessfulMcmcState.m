function saveUnsuccessfulMcmcState(storageDir, modelName, thinUsed)

if ~isfolder(storageDir)
  mkdir(storageDir);
end
lastUnsuccessfulThin = thinUsed;
save(fullfile(storageDir, sprintf('%s_mcmcState.mat', modelName)), ...
  'lastUnsuccessfulThin');
legacy = fullfile(storageDir, sprintf('%s_thinState.mat', modelName));
if isfile(legacy)
  delete(legacy);
end
end
