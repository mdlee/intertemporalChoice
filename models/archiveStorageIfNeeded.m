function archiveStorageIfNeeded(storageDir, modelName, dataName, engine, archiveTag)
%ARCHIVESTORAGEIFNEEDED  Copy canonical .mat to tagged path once (no overwrite of archive).
canonical = storageMatPath(storageDir, modelName, dataName, engine, '');
archived = storageMatPath(storageDir, modelName, dataName, engine, archiveTag);
if isfile(canonical) && ~isfile(archived)
  copyfile(canonical, archived);
  fprintf('Archived %s → %s\n', canonical, archived);
elseif isfile(archived)
  fprintf('Archive already present: %s\n', archived);
elseif ~isfile(canonical)
  fprintf('No canonical to archive: %s\n', canonical);
end
end
