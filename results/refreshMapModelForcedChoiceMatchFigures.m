function paths = refreshMapModelForcedChoiceMatchFigures()
%REFRESHMAPMODELFORCEDCHOICEMATCHFIGURES  Rebuild repo eps/png.
%
%   PATHS = refreshMapModelForcedChoiceMatchFigures()
%   Runs drawFiguresEntrop for analysis mapModelForcedChoiceMatch only,
%   writing:
%     results/descriptiveAdequacy/mapModelForcedChoiceMatch.eps
%     results/descriptiveAdequacy/mapModelForcedChoiceMatch.png
%   Safe to call from a long-running MCMC script (keeps caller workspace).
%   Half/double summaries rebuild automatically when storage .mat mtimes change.

resultsDir = fileparts(mfilename('fullpath'));
figuresDir = fullfile(resultsDir, 'descriptiveAdequacy');
epsPath = fullfile(figuresDir, 'mapModelForcedChoiceMatch.eps');
pngPath = fullfile(figuresDir, 'mapModelForcedChoiceMatch.png');

here = pwd;
cleanupPwd = onCleanup(@() cd(here));
cd(resultsDir);

analysisListOverride = {'mapModelForcedChoiceMatch'};
drawFiguresEntropKeepWorkspace = true;
run(fullfile(resultsDir, 'drawFiguresEntrop.m'));

addpath(fullfile(resultsDir, '..', 'models'));
addpath(fullfile(resultsDir, '..', 'general'));

if ~isfile(epsPath) || ~isfile(pngPath)
  error('refreshMapModelForcedChoiceMatchFigures:missingOutput', ...
    'Expected figure files were not written:\n  %s\n  %s', epsPath, pngPath);
end

paths = struct('eps', epsPath, 'png', pngPath, 'figuresDir', figuresDir);
end
