function fig = plotParticipantParameterCIs(chains, paramNames, varargin)
%PLOTPARTICIPANTPARAMETERCIS  5x5 panels: mean and 95% CI per parameter, per participant.
%   Each parameter uses one x-axis range and tick set shared across participants
%   (span = widest 95% CI over participants, plus padding).

p = inputParser;
addParameter(p, 'nRows', 5, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'nCols', 5, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'quantiles', [0.025 0.975], @(x) isnumeric(x) && numel(x) == 2);
addParameter(p, 'sgtitle', '', @ischar);
addParameter(p, 'savePath', '', @ischar);
addParameter(p, 'paramLabels', {}, @(x) iscell(x) || isstring(x));
addParameter(p, 'xPadFrac', 0.08, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});

paramNames = cellstr(paramNames(:));
nPar = numel(paramNames);
if isempty(p.Results.paramLabels)
  labels = paramNames;
else
  labels = cellstr(p.Results.paramLabels(:));
end

mats = cell(nPar, 1);
for k = 1:nPar
  mats{k} = codaIndexedMatrices(chains, paramNames{k});
end
nP = numel(mats{1});
nPanels = p.Results.nRows * p.Results.nCols;
q = p.Results.quantiles;
padFrac = p.Results.xPadFrac;

muAll = nan(nP, nPar);
qlAll = nan(nP, nPar);
qhAll = nan(nP, nPar);
for pp = 1:nP
  for k = 1:nPar
    s = mats{k}{pp}(:);
    muAll(pp, k) = mean(s);
    qlAll(pp, k) = quantile(s, q(1));
    qhAll(pp, k) = quantile(s, q(2));
  end
end

xLim = nan(nPar, 2);
xTicks = cell(nPar, 1);
for k = 1:nPar
  lo = inf;
  hi = -inf;
  for pp = 1:nP
    span = qhAll(pp, k) - qlAll(pp, k);
    if span <= 0 || ~isfinite(span)
      span = max(abs(muAll(pp, k)), 1e-6);
    end
    lo = min(lo, qlAll(pp, k) - padFrac * span);
    hi = max(hi, qhAll(pp, k) + padFrac * span);
  end
  if ~isfinite(lo) || ~isfinite(hi) || lo >= hi
    lo = -1;
    hi = 1;
  end
  xLim(k, :) = [lo hi];
  xTicks{k} = linspace(lo, hi, 5);
end

fig = figure('Color', 'w', 'Position', [60 60 1150 920]);
for pp = 1:nPanels
  if pp > nP
    subplot(p.Results.nRows, p.Results.nCols, pp);
    axis off;
    continue;
  end

  sp = subplot(p.Results.nRows, p.Results.nCols, pp);
  cellPos = get(sp, 'Position');
  delete(sp);

  rowGap = 0.02;
  rowH = (cellPos(4) - rowGap * (nPar - 1)) / nPar;
  for k = 1:nPar
    bottom = cellPos(2) + (k - 1) * (rowH + rowGap);
    ax = axes('Position', [cellPos(1), bottom, cellPos(3), rowH]);
    hold(ax, 'on');
    plot(ax, [qlAll(pp, k) qhAll(pp, k)], [0 0], '-', 'Color', [0.2 0.45 0.75], 'LineWidth', 2);
    plot(ax, muAll(pp, k), 0, 'o', 'MarkerFaceColor', [0.2 0.45 0.75], ...
      'MarkerEdgeColor', 'w', 'MarkerSize', 5);
    xlim(ax, xLim(k, :));
    xticks(ax, xTicks{k});
    xtickformat(ax, '%.2g');
    ylim(ax, [-1 1]);
    set(ax, 'YTick', [], 'TickDir', 'out', 'Box', 'off', 'FontSize', 6);
    if k == nPar
      title(ax, sprintf('p%d', pp), 'FontSize', 9, 'FontWeight', 'normal');
    end
    if mod(pp - 1, p.Results.nCols) == 0
      ylabel(ax, labels{k}, 'FontSize', 7);
    end
    hold(ax, 'off');
  end
end

if ~isempty(p.Results.sgtitle)
  topGap = 0.04;
  axAll = findobj(fig, 'Type', 'axes');
  for i = 1:numel(axAll)
    pos = axAll(i).Position;
    axAll(i).Position = [pos(1), max(pos(2) - topGap, 0.02), pos(3), pos(4)];
  end
  annotation(fig, 'textbox', [0 0.965 1 0.035], ...
    'String', p.Results.sgtitle, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontWeight', 'normal', ...
    'FontSize', 11, ...
    'FitBoxToText', 'off', ...
    'Interpreter', 'none');
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
