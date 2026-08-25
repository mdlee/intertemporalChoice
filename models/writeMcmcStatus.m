function writeMcmcStatus(statusFile, msg, varargin)
%WRITEMCMCSTATUS  Overwrite snapshot .status and append to running log.
%
%   writeMcmcStatus(statusFile, msg)
%   writeMcmcStatus(statusFile, msg, 'runningLog', path)
%
%   Snapshot file (statusFile): two lines — timestamp + message (overwrite).
%   Running log (default models/logs/mcmcRuns.status.log): append one line
%   "[timestamp] message" so progress can be tailed with:
%     tail -f models/logs/mcmcRuns.status.log

p = inputParser;
addParameter(p, 'runningLog', '', @(s) ischar(s) || isstring(s));
parse(p, varargin{:});

ts = datestr(now, 31);
msg = char(msg);

if ~isempty(statusFile)
  logDir = fileparts(statusFile);
  if ~isempty(logDir) && ~isfolder(logDir)
    mkdir(logDir);
  end
  fid = fopen(statusFile, 'w');
  if fid >= 0
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    fprintf(fid, '%s\n%s\n', ts, msg);
  end
end

runningLog = char(p.Results.runningLog);
if isempty(runningLog)
  modelsDir = fileparts(mfilename('fullpath'));
  runningLog = fullfile(modelsDir, 'logs', 'mcmcRuns.status.log');
end
logDir = fileparts(runningLog);
if ~isempty(logDir) && ~isfolder(logDir)
  mkdir(logDir);
end
fid2 = fopen(runningLog, 'a');
if fid2 >= 0
  cleaner2 = onCleanup(@() fclose(fid2)); %#ok<NASGU>
  fprintf(fid2, '[%s] %s\n', ts, msg);
end
drawnow;
try
  writeOvernightLiveCanvas();
catch
end
end
