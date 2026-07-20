function paths = refreshMapModelParameterRobustnessFigures()
%REFRESHMAPMODELPARAMETERROBUSTNESSFIGURES  Rebuild manuscript eps/png.
%
%   PATHS = refreshMapModelParameterRobustnessFigures()
%   Runs drawFiguresEntrop for analysis mapModelParameterRobustness only,
%   writing:
%     results/parameterInferences/mapModelParameterRobustness.eps
%     results/parameterInferences/mapModelParameterRobustness.png
%   Safe to call from a long-running MCMC script (keeps caller workspace).

resultsDir = fileparts(mfilename('fullpath'));
figuresDir = fullfile(resultsDir, 'parameterInferences');
epsPath = fullfile(figuresDir, 'mapModelParameterRobustness.eps');
pngPath = fullfile(figuresDir, 'mapModelParameterRobustness.png');

here = pwd;
cleanupPwd = onCleanup(@() cd(here));
cd(resultsDir);

analysisListOverride = {'mapModelParameterRobustness'};
drawFiguresEntropKeepWorkspace = true;
run(fullfile(resultsDir, 'drawFiguresEntrop.m'));

% drawFiguresEntrop may have left keep-workspace paths; ensure models/ is on path.
addpath(fullfile(resultsDir, '..', 'models'));
addpath(fullfile(resultsDir, '..', 'general'));

if ~isfile(epsPath) || ~isfile(pngPath)
  error('refreshMapModelParameterRobustnessFigures:missingOutput', ...
    'Expected figure files were not written:\n  %s\n  %s', epsPath, pngPath);
end

paths = struct('eps', epsPath, 'png', pngPath, 'figuresDir', figuresDir);
end
