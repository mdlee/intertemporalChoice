function fig = plotLatentMixtureZPosterior(chains, varargin)
%PLOTLATENTMIXTUREZPOSTERIOR  5x5 panels of P(model | data) per participant.
%
%   FIG = PLOTLATENTMIXTUREZPOSTERIOR(CHAINS) pools monitored z_1..z_n over
%   chains and bar-plots posterior model probabilities. One figure per call.
%
%   Name-value options:
%     nModels      — default 11
%     modelNames   — cell of short labels (default Ex..SS)
%     nRows, nCols — default 5, 5
%     sgtitle      — figure title
%     savePath     — if set, saveas(fig, savePath)

p = inputParser;
addParameter(p, 'nModels', 11, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'modelNames', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'nRows', 5, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'nCols', 5, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'sgtitle', '', @ischar);
addParameter(p, 'savePath', '', @ischar);
parse(p, varargin{:});
nModels = p.Results.nModels;

if isempty(p.Results.modelNames)
  modelNames = {'Ex'; 'Hc'; 'Hd'; 'PD'; 'DD'; 'Tr'; 'UT'; 'IT'; 'Gu'; 'LL'; 'SS'};
else
  modelNames = cellstr(p.Results.modelNames(:));
end

matrices = codaIndexedMatrices(chains, 'z');
nP = numel(matrices);
nPanels = p.Results.nRows * p.Results.nCols;

fig = figure('Color', 'w', 'Position', [80 80 1100 900]);
for pp = 1:nPanels
  subplot(p.Results.nRows, p.Results.nCols, pp);
  if pp > nP
    axis off;
    continue;
  end
  z = round(matrices{pp}(:));
  z(z < 1 | z > nModels) = [];
  if isempty(z)
    axis off;
    title(sprintf('p%d (no z)', pp), 'FontSize', 9);
    continue;
  end
  counts = histcounts(z, 0.5:(nModels + 0.5));
  probs = counts / sum(counts);
  bar(1:nModels, probs, 0.85, 'FaceColor', [0.35 0.55 0.85], 'EdgeColor', 'none');
  xlim([0.5 nModels + 0.5]);
  ylim([0 1]);
  set(gca, 'XTick', 1:nModels, 'XTickLabel', modelNames, 'XTickLabelRotation', 45, ...
    'FontSize', 7, 'TickDir', 'out', 'Box', 'off');
  title(sprintf('p%d', pp), 'FontSize', 9, 'FontWeight', 'normal');
  if mod(pp, p.Results.nCols) == 1
    ylabel('P(model|z)', 'FontSize', 8);
  end
end

if ~isempty(p.Results.sgtitle)
  sgtitle(p.Results.sgtitle, 'FontWeight', 'normal', 'FontSize', 11);
end

drawnow;

if ~isempty(p.Results.savePath)
  [saveDir, ~, ~] = fileparts(p.Results.savePath);
  if ~isempty(saveDir) && ~isfolder(saveDir)
    mkdir(saveDir);
  end
  saveas(fig, p.Results.savePath);
end

end
