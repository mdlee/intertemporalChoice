function writeOvernightLiveCanvas(canvasPath)
%WRITEOVERNIGHTLIVECANVAS  Live TSX canvas — Stan sequential in detail.
% Calls models/writeOvernightLiveCanvas.py so CmdStan CSV progress can be
% parsed while MATLAB is busy sampling.

if nargin < 1 || isempty(canvasPath)
  canvasPath = ['/home/mdlee/.cursor/projects/' ...
    'home-mdlee-Dropbox-GitHub-intertemporalChoice/canvases/' ...
    'overnight-mcmc-live.canvas.tsx'];
end

py = fullfile(fileparts(mfilename('fullpath')), 'writeOvernightLiveCanvas.py');
cmd = sprintf('python3 "%s" "%s"', py, canvasPath);
[st, out] = system(cmd);
if st ~= 0
  warning('writeOvernightLiveCanvas:python', '%s', strtrim(out));
end
end
