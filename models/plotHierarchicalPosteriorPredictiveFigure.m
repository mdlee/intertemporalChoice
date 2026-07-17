function fig = plotHierarchicalPosteriorPredictiveFigure(chains, data, d, modelName, varargin)
%PLOTHIERARCHICALPOSTERIORPREDICTIVEFIGURE  5x5 post. pred. count panels (one per p).
%
%   Same tile convention as plot_data_with_posterior_predictive_by_participant.m:
%   area proportional to P(replicated count = k); observed count highlighted.

p = inputParser;
addParameter(p, 'nRows', 5, @isnumeric);
addParameter(p, 'nCols', 5, @isnumeric);
addParameter(p, 'nPredSamples', 4000, @isnumeric);
addParameter(p, 'probDrawMin', [], @isnumeric);
addParameter(p, 'sgtitle', '', @ischar);
addParameter(p, 'savePath', '', @ischar);
parse(p, varargin{:});

nP = double(d.nParticipants);
nT = double(d.nTrials);
if isfield(d, 'nPairs')
  nPairs = double(d.nPairs);
else
  nPairs = double(max(d.pair(:)));
end

if isempty(p.Results.probDrawMin)
  probDrawMin = 1 / max(500, p.Results.nPredSamples);
else
  probDrawMin = p.Results.probDrawMin;
end

nTpMat = zeros(nP, nPairs);
obsCntMat = nan(nP, nPairs);
maxKByPar = zeros(1, nP);
for pp = 1:nP
  for pj = 1:nPairs
    m = d.pair(pp, :) == pj;
    if any(m)
      nTpMat(pp, pj) = sum(m);
      obsCntMat(pp, pj) = sum(d.LL(pp, m));
    end
  end
  maxKByPar(pp) = max(nTpMat(pp, :));
end

tileFace = [0.78 0.86 0.98];
tileEdge = [0.45 0.52 0.62];
obsFace = [0.92 0.48 0.22];
obsEdge = [0.45 0.2 0.08];
minObsSide = 0.36;

fig = figure('Color', 'w', 'Position', [40 40 1280 980], ...
  'Name', sprintf('%s: posterior predictive counts', modelName));

for pp = 1:(p.Results.nRows * p.Results.nCols)
  ax = subplot(p.Results.nRows, p.Results.nCols, pp);
  if pp > nP
    axis(ax, 'off');
    continue;
  end

  Th = hierarchicalExecutionThetaDraws(chains, data, modelName, pp);
  pmfPad = postPredictiveCountPmfFromTheta( ...
    Th, d.pair(pp, :), nPairs, nTpMat(pp, :), p.Results.nPredSamples);
  obsRow = obsCntMat(pp, :);
  maxK = maxKByPar(pp);

  if isempty(pmfPad) || all(isnan(pmfPad(:)), 'all')
    text(ax, 0.5, 0.5, 'No \theta', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    title(ax, sprintf('p%d', pp), 'FontWeight', 'normal');
    xlim(ax, [0.5 nPairs + 0.5]);
    ylim(ax, [-0.5 maxK + 0.5]);
    continue;
  end

  maxProb = max(pmfPad(:), [], 'omitnan');
  if ~(maxProb > 0)
    text(ax, 0.5, 0.5, 'No mass', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    title(ax, sprintf('p%d', pp), 'FontWeight', 'normal');
    continue;
  end

  scale = 0.92 / sqrt(maxProb);
  hold(ax, 'on');

  flat = pmfPad(:);
  [sortedProb, ord] = sort(flat, 'ascend', 'ComparisonMethod', 'real');
  [pjVec, kIdxVec] = ind2sub(size(pmfPad), ord);
  for z = 1:numel(sortedProb)
    prob = sortedProb(z);
    if isnan(prob) || prob < probDrawMin
      continue;
    end
    pj = pjVec(z);
    k = kIdxVec(z) - 1;
    nk = nTpMat(pp, pj);
    if nk <= 0 || k > nk || k < 0
      continue;
    end
    side = min(scale * sqrt(prob), 0.92);
    rectangle(ax, 'Position', [pj - side / 2, k - side / 2, side, side], ...
      'FaceColor', tileFace, 'EdgeColor', tileEdge, 'LineWidth', 0.65);
  end

  for pj = 1:nPairs
    nk = nTpMat(pp, pj);
    if nk <= 0
      continue;
    end
    ko = obsRow(pj);
    if ~isfinite(ko)
      continue;
    end
    ko = round(ko);
    if ko < 0 || ko > nk
      continue;
    end
    prob = pmfPad(pj, ko + 1);
    if isnan(prob) || prob < 0
      prob = 0;
    end
    rawSide = min(scale * sqrt(max(prob, eps)), 0.92);
    if prob < probDrawMin
      side = min(max(rawSide, minObsSide), 0.92);
    else
      side = rawSide;
    end
    rectangle(ax, 'Position', [pj - side / 2, ko - side / 2, side, side], ...
      'FaceColor', obsFace, 'EdgeColor', obsEdge, 'LineWidth', 1.5);
  end

  hold(ax, 'off');
  xlim(ax, [0.5 nPairs + 0.5]);
  ylim(ax, [-0.5 maxK + 0.5]);
  xlabel(ax, 'Problem', 'FontSize', 7);
  ylabel(ax, 'Count LL', 'FontSize', 7);
  title(ax, sprintf('p%d', pp), 'FontSize', 9, 'FontWeight', 'normal');
  set(ax, 'YTick', 0:maxK, 'XTick', 1:nPairs, 'FontSize', 6);
  grid(ax, 'on');
end

if ~isempty(p.Results.sgtitle)
  sgtitle(p.Results.sgtitle, 'FontWeight', 'normal', 'FontSize', 11);
else
  sgtitle(fig, sprintf('%s: posterior predictive counts (all participants)', modelName), ...
    'FontWeight', 'normal');
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
