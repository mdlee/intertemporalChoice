function [path, usedFallback] = resolveStorageMatPath(storageDir, modelName, dataName, engine, storageTag, verbose)
%RESOLVESTORAGEMATPATH  Tagged storage path, falling back to canonical if missing.
%
%   Models unchanged by an alternate run (e.g. cognitive fits under
%   storageTag='alphaPriorOld') were never archived as *__tag.mat. In that
%   case the canonical {model}_{data}_{engine}.mat is the correct prior-era
%   result. Set VERBOSE true to print a one-line notice when falling back.
%
%   See also: STORAGEMATPATH, ARCHIVESTORAGEIFNEEDED

if nargin < 5
  storageTag = '';
end
if nargin < 6 || isempty(verbose)
  verbose = false;
end
usedFallback = false;
path = storageMatPath(storageDir, modelName, dataName, engine, storageTag);
if isempty(char(storageTag)) || isfile(path)
  return;
end
canonical = storageMatPath(storageDir, modelName, dataName, engine, '');
if isfile(canonical)
  path = canonical;
  usedFallback = true;
  if verbose
    fprintf('  (no %s archive for %s — using canonical)\n', ...
      char(storageTag), modelName);
  end
end
end
