% Publication figures for hierarchical (entrop) intertemporal-choice analyses
%
% Toggle analyses in analysisList (comment in/out). Each case has graphics
% constants near the top for fine-tuning appearance.
% Global problem colors / display order are set once below (shared across figures).
%
% Output folders under results/:
%   basic/                experimentalDesign, dataCounts
%   latentMixture/        latentMixturePosteriors{,Half,Double},
%                         latentMixturePriorRobustness
%   descriptiveAdequacy/  mapModelPosteriorPredictive, allModelPosteriorPredictive,
%                         mapModelForcedChoiceMatch, modelAgreementTables
%   parameterInferences/  parameterInferences, allParameterInferences,
%                         mapModelParameterRobustness
%
% mapModelPosteriorPredictive caches a cutdown summary under models/storage/
% (mapModelPosteriorPredictive_summary_*.mat). Set regenerateMapPostPredSummaries
% inside that case to rebuild. allModelPosteriorPredictive caches
% allModelPosteriorPredictive_summary_*.mat (PP + forced-choice match for every
% model x participant). allParameterInferences writes one parameter-CI figure
% per model (all participants). modelAgreementTables writes APA LaTeX tables
% under manuscripts/ (forced-choice agreement + mean P(observed)).
% mapModelParameterRobustness and mapModelForcedChoiceMatch need half/double
% cognitive fits from models/runHierarchicalExecutionPriorRobustness.m.
%
% Run from results/:
%   drawFiguresEntrop
%
% From another script (e.g. incremental MCMC draft loop), set
%   analysisListOverride = {'mapModelParameterRobustness'};
%   % or {'mapModelForcedChoiceMatch'} / both
%   drawFiguresEntropKeepWorkspace = true;
% then run(this file) to refresh only those figures without clearing the caller.

keepWorkspace = exist('drawFiguresEntropKeepWorkspace', 'var') && ...
  logical(drawFiguresEntropKeepWorkspace);
overrideList = {};
if exist('analysisListOverride', 'var') && ~isempty(analysisListOverride)
  overrideList = analysisListOverride;
end
if ~keepWorkspace
  savedOverrideList = overrideList;
  clearvars -except savedOverrideList;
  close all;
  overrideList = savedOverrideList;
  clear savedOverrideList;
end

printFigures = true;

if ~isempty(overrideList)
  analysisList = overrideList;
else
  analysisList = {...
    %'experimentalDesign', ...
    %'dataCounts', ...
    %'latentMixturePosteriors', ...
    %'latentMixturePosteriorsHalf', ...
    %'latentMixturePosteriorsDouble', ...
    %'latentMixturePriorRobustness', ...
    %'mapModelPosteriorPredictive', ...
    %'allModelPosteriorPredictive', ...
    %'mapModelForcedChoiceMatch', ...
    %'modelAgreementTables', ...
    %'mapModelParameterRobustness', ...
    %'parameterInferences', ...
    'allParameterInferences', ...
    };
end
clear overrideList;

% ---- paths / data ----
resultsDir = fileparts(mfilename('fullpath'));
if isempty(resultsDir)
  resultsDir = pwd;
end
modelsDir = fullfile(resultsDir, '..', 'models');
generalDir = fullfile(resultsDir, '..', 'general');
storageDir = fullfile(modelsDir, 'storage');
addpath(modelsDir);
addpath(generalDir);

dataName = 'intertemporalChoice';
engine = 'jags';
mixtureName = 'latentMixtureHierarchicalPrecision_entrop';

% analysisName -> results subfolder
outputSubdirByAnalysis = containers.Map( ...
  {'experimentalDesign', 'dataCounts', ...
   'latentMixturePosteriors', 'latentMixturePosteriorsHalf', ...
   'latentMixturePosteriorsDouble', 'latentMixturePriorRobustness', ...
   'mapModelPosteriorPredictive', 'allModelPosteriorPredictive', ...
   'mapModelForcedChoiceMatch', 'modelAgreementTables', ...
   'mapModelParameterRobustness', 'parameterInferences', ...
   'allParameterInferences'}, ...
  {'basic', 'basic', ...
   'latentMixture', 'latentMixture', ...
   'latentMixture', 'latentMixture', ...
   'descriptiveAdequacy', 'descriptiveAdequacy', ...
   'descriptiveAdequacy', 'descriptiveAdequacy', ...
   'parameterInferences', 'parameterInferences', ...
   'parameterInferences'});

[data, d] = prepareIntertemporalChoiceData(dataName);
nP = double(d.nParticipants);
nT = double(d.nTrials);
if isfield(d, 'nPairs')
  nPairs = double(d.nPairs);
else
  nPairs = double(max(d.pair(:)));
end

% ---- shared graphics: problem colors (by pair type) ----
load pantoneColors pantone;
problemClr = {pantone.Freesia; ...
   pantone.CelosiaOrange; ...
   pantone.Paloma; ...
   pantone.PlacidBlue};

% Display order of original problem indices (1..nPairs).
% Affects numbering in experimentalDesign and x-axis order in data / post. pred.
% colors stay tied to each problem's pair type (unchanged by this permutation).
problemOrder = 1:nPairs;
% Example: problemOrder = [19:24, 1:18];

% Shared x-tick locations for problem axes (dataCounts, MAP post. pred., ...)
problemXtickLabels = [1 6 12 18 24];

if ~isequal(sort(problemOrder(:)'), 1:nPairs)
  error('problemOrder must be a permutation of 1:nPairs');
end
displayNumber = zeros(1, nPairs);
displayNumber(problemOrder) = 1:nPairs;

pairTypeByProblem = zeros(1, nPairs);
for pj = 1:nPairs
  m = find(d.pair == pj, 1);
  if isempty(m)
    error('No trials found for problem %d', pj);
  end
  pairTypeByProblem(pj) = d.pairType(m);
end

modelShort = {'Ex', 'Hc', 'Hd', 'PD', 'DD', 'Tr', 'UT', 'IT', 'Gu', 'LL', 'SS'};
modelLong = { ...
  'Exponential', 'Hyperbolic', 'Hyperboloid', 'Prop. differences', ...
  'Direct differences', 'Tradeoff', 'Unified tradeoff', 'ITCH', ...
  'Guess', 'LL response', 'SS response'};
nModels = numel(modelShort);

cognitiveStems = { ...
  'exponentialExecutionHierarchical_entrop', ...
  'hyperbolicExecutionHierarchical_entrop', ...
  'hyperboloidExecutionHierarchical_entrop', ...
  'proportionalDifferencesExecutionHierarchical_entrop', ...
  'directDifferencesExecutionHierarchical_entrop', ...
  'tradeoffExecutionHierarchical_entrop', ...
  'unifiedTradeoffExecutionHierarchical_entrop', ...
  'itchExecutionHierarchical_entrop'};
cognitiveJagsNames = cognitiveStems;
contaminantStems = {'guess', 'LL', 'SS'};

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
maxKGlobal = max(maxKByPar);

participantLabels = arrayfun(@(k) char('A' + k - 1), 1:nP, 'UniformOutput', false);

for analysisIdx = 1:numel(analysisList)
  analysisName = analysisList{analysisIdx};
  if ~isKey(outputSubdirByAnalysis, analysisName)
    error('No results subfolder mapped for analysis "%s"', analysisName);
  end
  figuresDir = fullfile(resultsDir, outputSubdirByAnalysis(analysisName));
  if ~isfolder(figuresDir)
    mkdir(figuresDir);
  end
  saveFigure = true;

  switch analysisName

    %% ================================================================
    case 'experimentalDesign'
      % Problem structure (from drawFigures_4), with display renumbering

      % graphics constants
      fontSize = 16;
      figPos = [0.2 0.1 0.45 0.85];
      height = 0.8;
      yLoc = [1, 2, 1, ...
        4, 5, ...
        7, ...
        1, 2, 1, ...
        4, 5, ...
        7, ...
        8, 9, 10, ...
        11, 12, 13, ...
        1, 2, 1, 3, 4, 5];
      labelColor = [0.5 0.5 0.5];

      F = figure; clf; hold on;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      subplot(2, 1, 1); cla; hold on;
      outcomesLabel = {'5150', '5400', '5450', '5600', '', '', '6050', '6200', '6350', '6500'};
      outcomes = [5150 5300 5450 5600 0 0 6050 6200 6350 6500];
      xTickLabel = {'1', '2', '3', '4', '', '', '7', '8', '9', '10'};
      set(gca, ...
        'xaxisloc', 'top', ...
        'ycolor', 'none', ...
        'xlim'       , [1 10]                , ...
        'xtick', 1:10 , ...
        'xticklabel', xTickLabel, ...
        'ylim'       , [0 12]                 , ...
        'ydir', 'rev', ...
        'tickdir', 'in', ...
        'ticklength', [0.01 0], ...
        'clipping', 'off', ...
        'fontsize'   , fontSize);
      moveAxis(gca, [1 1 1 1], [0.025 0 0 0]);
      H(1) = text(-0.25, 6, 'Large Outcomes', 'fontsize', fontSize+2, ...
        'vert', 'mid', 'hor', 'cen', 'rot', 90);
      H(2) =  text(0.5, 1.5, {'small', 'intervals'},  'fontsize', fontSize-2, ...
        'vert', 'mid', 'hor', 'cen','rot', 90);
      H(3) = text(0.5, 4.5, {'medium', 'intervals'},  'fontsize', fontSize-2, ...
        'vert', 'mid', 'hor', 'cen','rot', 90);
      H(4) = text(0.5, 9.5, {'large', 'intervals'},  'fontsize', fontSize-2, ...
        'vert', 'mid', 'hor', 'cen','rot', 90);

      H(5) =  text(-0.25, 2.5+17, 'Small Outcomes', 'fontsize', fontSize+2, ...
        'vert', 'mid', 'hor', 'cen', 'rot', 90);
      H(6) =  text(0.5, 2.5+17, {'large', 'intervals'},  'fontsize', fontSize-2, ...
        'vert', 'mid', 'hor', 'cen','rot', 90);

      set(H, 'color', labelColor);

      text(0.5, -1, 'outcome', ...
        'fontsize', fontSize, ...
        'hor', 'rig', 'ver', 'bot');
      text(0.5, -0.2, 'time', ...
        'fontsize', fontSize, ...
        'hor', 'rig', 'ver', 'bot');
      text(0.5, 15.6, 'outcome', ...
        'fontsize', fontSize, ...
        'hor', 'rig', 'ver', 'bot');
      text(0.5,16.4, 'time', ...
        'fontsize', fontSize, ...
        'hor', 'rig', 'ver', 'bot');
      for i = 1:10
        text(i, -1, outcomesLabel{i}, ...
          'fontsize', fontSize, ...
          'hor', 'cen', 'ver', 'bot');
      end

      for pIdx = 1:18
        match = find(d.pair == pIdx);
        rL = min(d.rewardA(match));
        rR = max(d.rewardA(match));
        rLidx = find(rL == outcomes);
        rRidx = find(rR == outcomes);
        pairType = pairTypeByProblem(pIdx);
        rectangle('position', [rLidx yLoc(pIdx)-height/2 rRidx-rLidx height], ...
          'curvature', [0 0], ...
          'facecolor', problemClr{pairType}, ...
          'edgecolor', 'none');
        text(rLidx+(rRidx-rLidx)/2, yLoc(pIdx), sprintf('%d', displayNumber(pIdx)), ...
          'fontsize', fontSize, ...
          'hor', 'cen', 'ver', 'mid');
      end

      subplot(2, 1, 2); cla; hold on;
      outcomesLabel = {'1150', '1250', '1350', '1450', '', '', '6050', '6200', '6350', '6500'};
      outcomes = [1150 1250 1350 1450];
      xTickLabel = {'1', '4', '7', '10'};
      set(gca, ...
        'xaxisloc', 'top', ...
        'ycolor', 'none', ...
        'xlim'       , [1 4]                , ...
        'xtick', 1:10 , ...
        'xticklabel', xTickLabel, ...
        'ylim'       , [0 12]                 , ...
        'ydir', 'rev', ...
        'tickdir', 'in', ...
        'ticklength', [0.01 0], ...
        'clipping', 'off', ...
        'fontsize'   , fontSize);
      moveAxis(gca, [1 1 1 1], [0.025 0 0 0]);

      for i = 1:10
        text(i, -1, outcomesLabel{i}, ...
          'fontsize', fontSize, ...
          'hor', 'cen', 'ver', 'bot');
      end
      for pIdx = 19:nPairs
        match = find(d.pair == pIdx);
        rL = min(d.rewardA(match));
        rR = max(d.rewardA(match));
        rLidx = find(rL == outcomes);
        rRidx = find(rR == outcomes);
        pairType = pairTypeByProblem(pIdx);
        rectangle('position', [rLidx yLoc(pIdx)-height/2 rRidx-rLidx height], ...
          'curvature', [0 0], ...
          'facecolor', problemClr{pairType}, ...
          'edgecolor', 'none');
        text(rLidx+(rRidx-rLidx)/2, yLoc(pIdx), sprintf('%d', displayNumber(pIdx)), ...
          'fontsize', fontSize, ...
          'hor', 'cen', 'ver', 'mid');
      end

    %% ================================================================
    case 'dataCounts'
      % 5x5: observed LL counts only (Kruschke tile layout)

      % graphics constants
      fontSize = 14;
      titleFontSize = 14;
      titleFontName = 'Helvetica';
      titleFontWeight = 'normal';
      labelFontSize = 16;
      nRows = 5;
      nCols = 5;
      figPos = [0.2 0.2 0.45 0.6];
      % One tile per (problem, count); axis units are integers so max side is 1.
      obsTileInset = 0.01; % clearance from neighboring tiles (per side)
      obsSide = 1 - 2 * obsTileInset;
      obsEdgeDarken = 0.35;
      obsLineWidth = 1.25;
      tickLength = 0.025;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 1 0.8];
      moveAxisShift = [0 0 0 0];
      showTitles = true;
      showXLabel = true;
      showYLabel = true;
      xLabelStr = 'Problem';
      yLabelStr = 'Count Larger-Later';
      % Outer label frame: [] = auto from content panels (ignores Raxes copies).
      % Or set explicitly, e.g. [0.08 0.06 0.86 0.88].
      supAxesPos = [];
      supAxesBuf = 0.05;
      supYLabelCloser = 0.025; % shift y superlabel right (toward panels)
      supXLabelCloser = 0.025; % shift x superlabel up (toward panels)
      showXTickLabels = true;
      showYTickLabels = true;
      tickLabelsOuterOnly = true;
      sharedYLim = true;
      showGrid = false;
      showSgtitle = false;
      sgtitleStr = 'Observed LL counts by problem';

      F = figure; clf;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      contentPos = nan(nRows * nCols, 4);
      for pp = 1:(nRows * nCols)
        ax = subplot(nRows, nCols, pp);
        if pp > nP
          axis(ax, 'off');
          continue;
        end

        hold(ax, 'on');
        obsRow = obsCntMat(pp, :);
        for dispIdx = 1:nPairs
          pj = problemOrder(dispIdx);
          ko = obsRow(pj);
          if ~isfinite(ko)
            continue;
          end
          ko = round(ko);
          face = problemClr{pairTypeByProblem(pj)};
          edge = max(0, face * (1 - obsEdgeDarken));
          rectangle(ax, 'Position', [dispIdx - obsSide / 2, ko - obsSide / 2, obsSide, obsSide], ...
            'FaceColor', face, 'EdgeColor', edge, 'LineWidth', obsLineWidth);
        end
        hold(ax, 'off');

        if sharedYLim
          yHi = maxKGlobal;
        else
          yHi = maxKByPar(pp);
        end
        xlim(ax, [1 nPairs]);
        ylim(ax, [0 yHi]);
        set(ax, ...
          'TickDir', 'out', ...
          'Box', 'off', ...
          'FontSize', fontSize, ...
          'XTick', problemXtickLabels, ...
          'xticklabelrot', 0, ...
          'YTick', [0 yHi], ...
          'TickLength', [tickLength 0], ...
          'clipping', 'off');
        onBottom = pp > (nRows - 1) * nCols;
        onLeft = mod(pp, nCols) == 1;
        showXHere = showXTickLabels && (~tickLabelsOuterOnly || onBottom);
        showYHere = showYTickLabels && (~tickLabelsOuterOnly || onLeft);
        if ~showXHere
          set(ax, 'XTickLabel', []);
        end
        if ~showYHere
          set(ax, 'YTickLabel', []);
        end
        if showGrid
          grid(ax, 'on');
        end
        moveAxis(gca, moveAxisScale, moveAxisShift);
        contentPos(pp, :) = get(ax, 'Position');
        [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
        set([axX axY], 'Tag', 'RaxesCopy');
        % Raxes copyobj keeps titles on the shifted axes; keep label only on content ax
        title(axX, '');
        title(axY, '');
        if showTitles
          ht = title(ax, participantLabels{pp}, ...
            'FontName', titleFontName, ...
            'FontSize', titleFontSize, ...
            'FontWeight', titleFontWeight, ...
            'HorizontalAlignment', 'left');
          xl = xlim(ax);
          tp = get(ht, 'Position');
          set(ht, 'Position', [xl(1), tp(2), tp(3)]);
        else
          title(ax, '');
        end
      end

      if isempty(supAxesPos)
        keep = all(isfinite(contentPos), 2);
        pos = contentPos(keep, :);
        leftMin = min(pos(:, 1));
        bottomMin = min(pos(:, 2));
        leftMax = max(pos(:, 1) + pos(:, 3));
        bottomMax = max(pos(:, 2) + pos(:, 4));
        % Expand for outward Raxes shifts when placing outer labels
        supAxesPos = [ ...
          leftMin - supAxesBuf - raxesYShift, ...
          bottomMin - supAxesBuf - raxesXShift, ...
          (leftMax - leftMin) + 2 * supAxesBuf + raxesYShift, ...
          (bottomMax - bottomMin) + 2 * supAxesBuf + raxesXShift];
      end

      if showXLabel
        [~, Hx] = suplabel(xLabelStr, 'x', supAxesPos);
        hxPos = get(Hx, 'Position');
        hxPos(2) = hxPos(2) + supXLabelCloser;
        set(Hx, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
          'Position', hxPos);
      end
      if showYLabel
        [~, Hy] = suplabel(yLabelStr, 'y', supAxesPos);
        hyPos = get(Hy, 'Position');
        hyPos(1) = hyPos(1) + supYLabelCloser;
        set(Hy, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
          'Position', hyPos);
      end

      if showSgtitle
        sgtitle(sgtitleStr, 'FontWeight', 'normal', 'FontSize', titleFontSize + 2);
      end

    %% ================================================================
    case {'latentMixturePosteriors', 'latentMixturePosteriorsHalf', 'latentMixturePosteriorsDouble'}
      % 5x5: P(model | data) from hierarchical latent mixture (entrop)
      switch analysisName
        case 'latentMixturePosteriorsHalf'
          mixStem = 'latentMixtureHierarchicalMuPrecHalf_entrop';
          sgtitleStr = 'Latent mixture posteriors (half priors)';
        case 'latentMixturePosteriorsDouble'
          mixStem = 'latentMixtureHierarchicalMuPrecDouble_entrop';
          sgtitleStr = 'Latent mixture posteriors (double priors)';
        otherwise
          mixStem = mixtureName; % latentMixtureHierarchicalPrecision_entrop
          sgtitleStr = 'Latent mixture posterior model probabilities';
      end

      % graphics constants
      fontSize = 10;
      titleFontSize = 14;
      titleFontName = 'Helvetica';
      titleFontWeight = 'normal';
      labelFontSize = 16;
      nRows = 5;
      nCols = 5;
      figPos = [0.2 0.2 0.45 0.6];
      barFace = [0.35 0.55 0.85];
      barEdge = 'none';
      barWidth = 0.85;
      tickLength = 0.025;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 1 0.8];
      moveAxisShift = [0 0 0 0];
      showTitles = true;
      showXLabel = true;
      showYLabel = true;
      xLabelStr = 'Model';
      yLabelStr = 'Posterior Model Probability';
      % Outer label frame: [] = auto from content panels (ignores Raxes copies).
      % Or set explicitly, e.g. [0.08 0.06 0.86 0.88].
      supAxesPos = [];
      supAxesBuf = 0.05;
      supYLabelCloser = 0.025; % shift y superlabel right (toward panels)
      supXLabelCloser = 0.025; % shift x superlabel up (toward panels)
      xtickLabels = 1:nModels;
      xTickLabelRotation = 90;
      showXTickLabels = true;
      showYTickLabels = true;
      tickLabelsOuterOnly = true;
      showGrid = false;
      showSgtitle = false;
      yLim = [0 1];
      yTicks = [0 1];

      mixPath = fullfile(storageDir, sprintf('%s_%s_%s.mat', mixStem, dataName, engine));
      if ~isfile(mixPath)
        error('Missing mixture chains: %s', mixPath);
      end
      fprintf('Loading %s\n', mixPath);
      mixS = load(mixPath, 'chains');
      zMats = codaIndexedMatrices(mixS.chains, 'z');

      F = figure; clf;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      contentPos = nan(nRows * nCols, 4);
      for pp = 1:(nRows * nCols)
        ax = subplot(nRows, nCols, pp);
        if pp > nP || pp > numel(zMats)
          axis(ax, 'off');
          continue;
        end

        z = round(zMats{pp}(:));
        z = z(isfinite(z) & z >= 1 & z <= nModels);
        if isempty(z)
          axis(ax, 'off');
          if showTitles
            ht = title(ax, sprintf('%s (no z)', participantLabels{pp}), ...
              'FontName', titleFontName, ...
              'FontSize', titleFontSize, ...
              'FontWeight', titleFontWeight, ...
              'HorizontalAlignment', 'left');
            xl = xlim(ax);
            tp = get(ht, 'Position');
            set(ht, 'Position', [xl(1), tp(2), tp(3)]);
          end
          continue;
        end

        counts = histcounts(z, 0.5:(nModels + 0.5));
        probs = counts / sum(counts);
        bar(ax, 1:nModels, probs, barWidth, ...
          'FaceColor', barFace, 'EdgeColor', barEdge);
        xlim(ax, [1 nModels]);
        ylim(ax, yLim);
        set(ax, ...
          'TickDir', 'out', ...
          'Box', 'off', ...
          'FontSize', fontSize, ...
          'XTick', xtickLabels, ...
          'xticklabelrot', xTickLabelRotation, ...
          'YTick', yTicks, ...
          'TickLength', [tickLength 0], ...
          'clipping', 'off');
        onBottom = pp > (nRows - 1) * nCols;
        onLeft = mod(pp, nCols) == 1;
        showXHere = showXTickLabels && (~tickLabelsOuterOnly || onBottom);
        showYHere = showYTickLabels && (~tickLabelsOuterOnly || onLeft);
        if showXHere
          set(ax, 'XTickLabel', lower(modelShort));
        else
          set(ax, 'XTickLabel', []);
        end
        if ~showYHere
          set(ax, 'YTickLabel', []);
        end
        if showGrid
          grid(ax, 'on');
        end
        moveAxis(gca, moveAxisScale, moveAxisShift);
        contentPos(pp, :) = get(ax, 'Position');
        [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
        set([axX axY], 'Tag', 'RaxesCopy');
        title(axX, '');
        title(axY, '');
        if showTitles
          ht = title(ax, participantLabels{pp}, ...
            'FontName', titleFontName, ...
            'FontSize', titleFontSize, ...
            'FontWeight', titleFontWeight, ...
            'HorizontalAlignment', 'left');
          xl = xlim(ax);
          tp = get(ht, 'Position');
          set(ht, 'Position', [xl(1), tp(2), tp(3)]);
        else
          title(ax, '');
        end
      end

      if isempty(supAxesPos)
        keep = all(isfinite(contentPos), 2);
        pos = contentPos(keep, :);
        leftMin = min(pos(:, 1));
        bottomMin = min(pos(:, 2));
        leftMax = max(pos(:, 1) + pos(:, 3));
        bottomMax = max(pos(:, 2) + pos(:, 4));
        supAxesPos = [ ...
          leftMin - supAxesBuf - raxesYShift, ...
          bottomMin - supAxesBuf - raxesXShift, ...
          (leftMax - leftMin) + 2 * supAxesBuf + raxesYShift, ...
          (bottomMax - bottomMin) + 2 * supAxesBuf + raxesXShift];
      end

      if showXLabel
        [~, Hx] = suplabel(xLabelStr, 'x', supAxesPos);
        hxPos = get(Hx, 'Position');
        hxPos(2) = hxPos(2) + supXLabelCloser;
        set(Hx, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
          'Position', hxPos);
      end
      if showYLabel
        [~, Hy] = suplabel(yLabelStr, 'y', supAxesPos);
        hyPos = get(Hy, 'Position');
        hyPos(1) = hyPos(1) + supYLabelCloser;
        set(Hy, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
          'Position', hyPos);
      end

      if showSgtitle
        sgtitle(sgtitleStr, 'FontWeight', 'normal', 'FontSize', titleFontSize + 2);
      end

    %% ================================================================
    case 'latentMixturePriorRobustness'
      % Console MAP stability under half / double priors, plus scatter of
      % posterior model probability (orig vs half / orig vs double) for all
      % model families (cognitive + contaminant), with Tr and UT combined.

      mixStems = { ...
        mixtureName, ...
        'latentMixtureHierarchicalMuPrecHalf_entrop', ...
        'latentMixtureHierarchicalMuPrecDouble_entrop'};
      mixLabels = {'original', 'half', 'double'};
      Pcell = cell(1, 3);
      for k = 1:3
        mixPath = fullfile(storageDir, sprintf('%s_%s_%s.mat', mixStems{k}, dataName, engine));
        if ~isfile(mixPath)
          error('Missing mixture chains (%s): %s', mixLabels{k}, mixPath);
        end
        fprintf('Loading %s\n', mixPath);
        mixS = load(mixPath, 'chains');
        Pcell{k} = posteriorModelProbsFromZ(mixS.chains, nP, nModels);
        clear mixS;
      end
      Porig = Pcell{1};
      Phalf = Pcell{2};
      Pdouble = Pcell{3};

      [~, mapModelOrig] = max(Porig, [], 1);
      [~, mapModelHalf] = max(Phalf, [], 1);
      [~, mapModelDouble] = max(Pdouble, [], 1);

      sameHalf = sum(mapModelOrig == mapModelHalf);
      sameDouble = sum(mapModelOrig == mapModelDouble);
      fprintf('\nMAP model stability under altered priors (n = %d participants):\n', nP);
      fprintf('  Original vs half:   %d / %d unchanged (%.0f%%)\n', ...
        sameHalf, nP, 100 * sameHalf / nP);
      fprintf('  Original vs double: %d / %d unchanged (%.0f%%)\n', ...
        sameDouble, nP, 100 * sameDouble / nP);
      changedHalf = find(mapModelOrig ~= mapModelHalf);
      if ~isempty(changedHalf)
        fprintf('  Changed under half:');
        for ii = 1:numel(changedHalf)
          pp = changedHalf(ii);
          fprintf(' %s(%s→%s)', participantLabels{pp}, ...
            lower(modelShort{mapModelOrig(pp)}), lower(modelShort{mapModelHalf(pp)}));
        end
        fprintf('\n');
      end
      changedDouble = find(mapModelOrig ~= mapModelDouble);
      if ~isempty(changedDouble)
        fprintf('  Changed under double:');
        for ii = 1:numel(changedDouble)
          pp = changedDouble(ii);
          fprintf(' %s(%s→%s)', participantLabels{pp}, ...
            lower(modelShort{mapModelOrig(pp)}), lower(modelShort{mapModelDouble(pp)}));
        end
        fprintf('\n');
      end

      % Treat Tr (6) and UT (7) as the same model-use inference
      mapFamilyOrig = mapModelToUseFamily(mapModelOrig);
      mapFamilyHalf = mapModelToUseFamily(mapModelHalf);
      mapFamilyDouble = mapModelToUseFamily(mapModelDouble);
      sameHalfFam = sum(mapFamilyOrig == mapFamilyHalf);
      sameDoubleFam = sum(mapFamilyOrig == mapFamilyDouble);
      fprintf('\nMAP stability treating tr/ut as the same model use (n = %d):\n', nP);
      fprintf('  Original vs half:   %d / %d unchanged (%.0f%%)\n', ...
        sameHalfFam, nP, 100 * sameHalfFam / nP);
      fprintf('  Original vs double: %d / %d unchanged (%.0f%%)\n', ...
        sameDoubleFam, nP, 100 * sameDoubleFam / nP);
      changedHalfFam = find(mapFamilyOrig ~= mapFamilyHalf);
      if ~isempty(changedHalfFam)
        fprintf('  Changed under half:');
        for ii = 1:numel(changedHalfFam)
          pp = changedHalfFam(ii);
          fprintf(' %s(%s→%s)', participantLabels{pp}, ...
            lower(modelShort{mapModelOrig(pp)}), lower(modelShort{mapModelHalf(pp)}));
        end
        fprintf('\n');
      end
      changedDoubleFam = find(mapFamilyOrig ~= mapFamilyDouble);
      if ~isempty(changedDoubleFam)
        fprintf('  Changed under double:');
        for ii = 1:numel(changedDoubleFam)
          pp = changedDoubleFam(ii);
          fprintf(' %s(%s→%s)', participantLabels{pp}, ...
            lower(modelShort{mapModelOrig(pp)}), lower(modelShort{mapModelDouble(pp)}));
        end
        fprintf('\n');
      end
      fprintf('\n');

      % Combined Tr+UT into one family; scatter P_orig vs P_half/double for
      % every family (all cognitive + contaminant), every participant.
      % Label points farther than dCritical from the diagonal with
      % MAP_orig->MAP_alt (combined-family two-letter codes).
      PfamOrig = combinedFamilyPosteriors(Porig);
      PfamHalf = combinedFamilyPosteriors(Phalf);
      PfamDouble = combinedFamilyPosteriors(Pdouble);
      nFam = size(PfamOrig, 1);
      familyShort = {'Ex', 'Hc', 'Hd', 'PD', 'DD', 'Tr', 'IT', 'Gu', 'LL', 'SS'};
      [~, mapFamOrig] = max(PfamOrig, [], 1);
      [~, mapFamHalf] = max(PfamHalf, [], 1);
      [~, mapFamDouble] = max(PfamDouble, [], 1);
      mapFamAlt = {mapFamHalf, mapFamDouble};
      fprintf('Scatter: %d participants × %d families (Tr+UT combined) = %d points per panel\n\n', ...
        nP, nFam, nP * nFam);

      fontSize = 14;
      labelFontSize = 16;
      titleFontSize = 14;
      pointLabelFontSize = 8;
      figPos = [0.15 0.25 0.7 0.45];
      markerSize = 36;
      markerFace = [0.45 0.55 0.75];
      markerEdge = 0.25 * [1 1 1];
      identityClr = 0.55 * [1 1 1];
      identityLineWidth = 1.0;
      tickLength = 0.02;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 0.9 0.9];
      moveAxisShift = [0 0.025 0 0];
      xLabelStr = 'Original prior';
      yLabelHalf = 'Half prior';
      yLabelDouble = 'Double prior';
      showPanelTitles = true;
      panelTitleHalf = 'Half';
      panelTitleDouble = 'Double';
      dCritical = 0.3; % |P_orig - P_alt| threshold for labeling off-diagonal points
      pointLabelOffset = 0.02;

      altP = {PfamHalf, PfamDouble};
      yLabs = {yLabelHalf, yLabelDouble};
      panTitles = {panelTitleHalf, panelTitleDouble};

      F = figure; clf;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      for pan = 1:2
        ax = subplot(1, 2, pan);
        hold(ax, 'on');
        plot(ax, [0 1], [0 1], '-', 'Color', identityClr, 'LineWidth', identityLineWidth);

        for fam = 1:nFam
          for pp = 1:nP
            xo = PfamOrig(fam, pp);
            ya = altP{pan}(fam, pp);
            if ~(isfinite(xo) && isfinite(ya))
              continue;
            end
            scatter(ax, xo, ya, markerSize, ...
              'MarkerFaceColor', markerFace, ...
              'MarkerEdgeColor', markerEdge, ...
              'LineWidth', 0.6);
          end
        end

        % One label per participant at the farthest-from-diagonal family point
        for pp = 1:nP
          dFam = abs(PfamOrig(:, pp) - altP{pan}(:, pp));
          [dMax, famStar] = max(dFam);
          if ~(isfinite(dMax) && dMax > dCritical)
            continue;
          end
          xo = PfamOrig(famStar, pp);
          ya = altP{pan}(famStar, pp);
          lbl = sprintf('%s->%s', ...
            lower(familyShort{mapFamOrig(pp)}), ...
            lower(familyShort{mapFamAlt{pan}(pp)}));
          if ya >= xo
            tx = xo - pointLabelOffset;
            ty = ya + pointLabelOffset;
            hAlign = 'right';
            vAlign = 'bottom';
          else
            tx = xo + pointLabelOffset;
            ty = ya - pointLabelOffset;
            hAlign = 'left';
            vAlign = 'top';
          end
          text(ax, tx, ty, lbl, ...
            'FontSize', pointLabelFontSize, ...
            'FontWeight', 'normal', ...
            'HorizontalAlignment', hAlign, ...
            'VerticalAlignment', vAlign, ...
            'Clipping', 'off', ...
            'Interpreter', 'none');
        end

        xlim(ax, [0 1]);
        ylim(ax, [0 1]);
        axis(ax, 'square');
        set(ax, ...
          'TickDir', 'out', ...
          'Box', 'off', ...
          'FontSize', fontSize, ...
          'XTick', 0:0.25:1, ...
          'YTick', 0:0.25:1, ...
          'TickLength', [tickLength 0], ...
          'clipping', 'off');
        xlabel(ax, xLabelStr, 'FontSize', labelFontSize);
        ylabel(ax, yLabs{pan}, 'FontSize', labelFontSize);
        if showPanelTitles
          title(ax, panTitles{pan}, 'FontSize', titleFontSize, 'FontWeight', 'normal');
        end
        moveAxis(gca, moveAxisScale, moveAxisShift);
        [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
        set([axX axY], 'Tag', 'RaxesCopy');
        title(axX, '');
        title(axY, '');
        if showPanelTitles
          title(ax, panTitles{pan}, 'FontSize', titleFontSize, 'FontWeight', 'normal');
        end
      end

      saveFigure = true;

    %% ================================================================
    case 'mapModelPosteriorPredictive'
      % 5x5: post. pred. under MAP model from latent mixture

      % graphics constants
      fontSize = 14;
      titleFontSize = 14;
      titleFontName = 'Helvetica';
      titleFontWeight = 'normal';
      labelFontSize = 16;
      nRows = 5;
      nCols = 5;
      figPos = [0.2 0.2 0.45 0.6];
      nPredSamples = 4000;
      probDrawMin = 1 / max(500, nPredSamples);
      tileLineWidth = 0.65;
      obsEdgeDarken = 0.35;
      obsLineWidth = 1.5;
      minObsSide = 0.36;
      predFaceClr = [1 1 1];
      predEdgeClr = 0.55 * [1 1 1];
      tickLength = 0.025;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 1 0.8];
      moveAxisShift = [0 0 0 0];
      showTitles = true;
      showModelInTitle = true;
      showMatchInTitle = true; % forced-choice agreement % (not posterior model prob)
      choiceThreshold = 0.5;
      modelTitleFontSize = 14;
      modelTitleFontName = 'Helvetica';
      modelTitleFontWeight = 'normal';
      titleYNorm = 1.06; % normalized axes units (1 = top of plot box)
      showXLabel = true;
      showYLabel = true;
      xLabelStr = 'Problem';
      yLabelStr = 'Count LL';
      % Outer label frame: [] = auto from content panels (ignores Raxes copies).
      % Or set explicitly, e.g. [0.08 0.06 0.86 0.88].
      supAxesPos = [];
      supAxesBuf = 0.05;
      supYLabelCloser = 0.025; % shift y superlabel right (toward panels)
      supXLabelCloser = 0.025; % shift x superlabel up (toward panels)
      showXTickLabels = true;
      showYTickLabels = true;
      tickLabelsOuterOnly = true;
      sharedYLim = true;
      showGrid = false;
      showSgtitle = false;
      sgtitleStr = 'MAP-model posterior predictive counts';
      rngSeed = 42;
      % Cutdown summary: mapModel / mapProb / pmfPad only (full chain mats kept).
      % Set true after refitting, or leave false to auto-rebuild when sources change.
      regenerateMapPostPredSummaries = false;
      regenerateMapForcedChoiceSummaries = false;

      mixPath = fullfile(storageDir, sprintf('%s_%s_%s.mat', mixtureName, dataName, engine));
      if ~isfile(mixPath)
        error('Missing mixture chains: %s', mixPath);
      end
      summaryPath = fullfile(storageDir, ...
        sprintf('mapModelPosteriorPredictive_summary_%s_%s.mat', dataName, engine));

      [mapModel, ~, pmfPadByParticipant, rebuilt] = loadOrBuildMapPostPredSummary( ...
        summaryPath, regenerateMapPostPredSummaries, ...
        mixPath, storageDir, dataName, engine, ...
        mixtureName, cognitiveStems, contaminantStems, cognitiveJagsNames, modelLong, ...
        data, d, nP, nT, nPairs, nModels, nTpMat, nPredSamples, rngSeed);
      if rebuilt
        fprintf('Wrote MAP post-pred summary %s\n', summaryPath);
      else
        fprintf('Loaded MAP post-pred summary %s\n', summaryPath);
      end

      matchSummaryPath = fullfile(storageDir, ...
        sprintf('mapModelForcedChoiceMatch_summary_%s_%s_%s.mat', ...
          mixtureName, dataName, engine));
      [~, ~, matchProp, ~, matchRebuilt] = loadOrBuildMapForcedChoiceMatchSummary( ...
        matchSummaryPath, regenerateMapForcedChoiceSummaries, ...
        mixPath, storageDir, dataName, engine, ...
        mixtureName, cognitiveStems, contaminantStems, cognitiveJagsNames, modelLong, ...
        data, d, nP, nT, nModels, choiceThreshold);
      if matchRebuilt
        fprintf('Wrote forced-choice match summary %s\n', matchSummaryPath);
      else
        fprintf('Loaded forced-choice match summary %s\n', matchSummaryPath);
      end

      F = figure; clf;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      contentPos = nan(nRows * nCols, 4);
      for pp = 1:(nRows * nCols)
        ax = subplot(nRows, nCols, pp);
        if pp > nP
          axis(ax, 'off');
          continue;
        end

        mi = mapModel(pp);
        if ~(isfinite(mi) && mi >= 1 && mi <= nModels) || isempty(pmfPadByParticipant{pp})
          text(ax, 0.5, 0.5, 'No fit', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', fontSize);
          if showTitles
            text(ax, 0, titleYNorm, participantLabels{pp}, ...
              'Units', 'normalized', ...
              'HorizontalAlignment', 'left', ...
              'VerticalAlignment', 'bottom', ...
              'FontName', titleFontName, ...
              'FontSize', titleFontSize, ...
              'FontWeight', titleFontWeight, ...
              'Interpreter', 'none', ...
              'Clipping', 'off');
          end
          continue;
        end

        pmfPad = pmfPadByParticipant{pp};

        if sharedYLim
          yHi = maxKGlobal;
        else
          yHi = maxKByPar(pp);
        end

        if isempty(pmfPad) || all(isnan(pmfPad(:)), 'all')
          text(ax, 0.5, 0.5, 'No \theta', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', fontSize);
        else
          maxProb = max(pmfPad(:), [], 'omitnan');
          if ~(maxProb > 0)
            text(ax, 0.5, 0.5, 'No mass', 'Units', 'normalized', ...
              'HorizontalAlignment', 'center', 'FontSize', fontSize);
          else
            scale = 0.92 / sqrt(maxProb);
            hold(ax, 'on');

            obsRow = obsCntMat(pp, :);
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
              dispIdx = displayNumber(pj);
              side = min(scale * sqrt(prob), 0.92);
              ko = obsRow(pj);
              isObs = isfinite(ko) && round(ko) == k;
              if isObs
                continue; % drawn below with problem color + thicker edge
              end
              rectangle(ax, 'Position', [dispIdx - side / 2, k - side / 2, side, side], ...
                'FaceColor', predFaceClr, 'EdgeColor', predEdgeClr, 'LineWidth', tileLineWidth);
            end

            for dispIdx = 1:nPairs
              pj = problemOrder(dispIdx);
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
              face = problemClr{pairTypeByProblem(pj)};
              edge = max(0, face * (1 - obsEdgeDarken));
              rectangle(ax, 'Position', [dispIdx - side / 2, ko - side / 2, side, side], ...
                'FaceColor', face, 'EdgeColor', edge, 'LineWidth', obsLineWidth);
            end
            hold(ax, 'off');
          end
        end

        xlim(ax, [1 nPairs]);
        ylim(ax, [0 yHi]);
        set(ax, ...
          'TickDir', 'out', ...
          'Box', 'off', ...
          'FontSize', fontSize, ...
          'XTick', problemXtickLabels, ...
          'xticklabelrot', 0, ...
          'YTick', [0 yHi], ...
          'TickLength', [tickLength 0], ...
          'clipping', 'off');
        onBottom = pp > (nRows - 1) * nCols;
        onLeft = mod(pp, nCols) == 1;
        showXHere = showXTickLabels && (~tickLabelsOuterOnly || onBottom);
        showYHere = showYTickLabels && (~tickLabelsOuterOnly || onLeft);
        if ~showXHere
          set(ax, 'XTickLabel', []);
        end
        if ~showYHere
          set(ax, 'YTickLabel', []);
        end
        if showGrid
          grid(ax, 'on');
        end
        moveAxis(gca, moveAxisScale, moveAxisShift);
        contentPos(pp, :) = get(ax, 'Position');
        [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
        set([axX axY], 'Tag', 'RaxesCopy');
        title(axX, '');
        title(axY, '');
        title(ax, '');
        if showTitles
          text(ax, 0, titleYNorm, participantLabels{pp}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', titleFontName, ...
            'FontSize', titleFontSize, ...
            'FontWeight', titleFontWeight, ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
          modelStr = '';
          matchPct = NaN;
          if showMatchInTitle && isfinite(matchProp(pp))
            matchPct = round(100 * matchProp(pp));
          end
          if showModelInTitle && showMatchInTitle && isfinite(matchPct)
            modelStr = sprintf('%s[%d%%]', lower(modelShort{mi}), matchPct);
          elseif showModelInTitle
            modelStr = lower(modelShort{mi});
          elseif showMatchInTitle && isfinite(matchPct)
            modelStr = sprintf('[%d%%]', matchPct);
          end
          if ~isempty(modelStr)
            text(ax, 1, titleYNorm, modelStr, ...
              'Units', 'normalized', ...
              'HorizontalAlignment', 'right', ...
              'VerticalAlignment', 'bottom', ...
              'FontName', modelTitleFontName, ...
              'FontSize', modelTitleFontSize, ...
              'FontWeight', modelTitleFontWeight, ...
              'Interpreter', 'none', ...
              'Clipping', 'off');
          end
        end
      end

      if isempty(supAxesPos)
        keep = all(isfinite(contentPos), 2);
        pos = contentPos(keep, :);
        leftMin = min(pos(:, 1));
        bottomMin = min(pos(:, 2));
        leftMax = max(pos(:, 1) + pos(:, 3));
        bottomMax = max(pos(:, 2) + pos(:, 4));
        supAxesPos = [ ...
          leftMin - supAxesBuf - raxesYShift, ...
          bottomMin - supAxesBuf - raxesXShift, ...
          (leftMax - leftMin) + 2 * supAxesBuf + raxesYShift, ...
          (bottomMax - bottomMin) + 2 * supAxesBuf + raxesXShift];
      end

      if showXLabel
        [~, Hx] = suplabel(xLabelStr, 'x', supAxesPos);
        hxPos = get(Hx, 'Position');
        hxPos(2) = hxPos(2) + supXLabelCloser;
        set(Hx, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
          'Position', hxPos);
      end
      if showYLabel
        [~, Hy] = suplabel(yLabelStr, 'y', supAxesPos);
        hyPos = get(Hy, 'Position');
        hyPos(1) = hyPos(1) + supYLabelCloser;
        set(Hy, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
          'Position', hyPos);
      end

      if showSgtitle
        sgtitle(sgtitleStr, 'FontWeight', 'normal', 'FontSize', titleFontSize + 2);
      end

    %% ================================================================
    case 'allModelPosteriorPredictive'
      % One 5x5 figure per model (8 cognitive + 3 contaminants): posterior
      % predictive counts for that model applied to all participants. Panel
      % labels show forced-choice agreement for that model.

      % graphics constants (match mapModelPosteriorPredictive)
      fontSize = 14;
      titleFontSize = 14;
      titleFontName = 'Helvetica';
      titleFontWeight = 'normal';
      labelFontSize = 16;
      nRows = 5;
      nCols = 5;
      figPos = [0.2 0.2 0.45 0.6];
      nPredSamples = 4000;
      probDrawMin = 1 / max(500, nPredSamples);
      tileLineWidth = 0.65;
      obsEdgeDarken = 0.35;
      obsLineWidth = 1.5;
      minObsSide = 0.36;
      predFaceClr = [1 1 1];
      predEdgeClr = 0.55 * [1 1 1];
      tickLength = 0.025;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 1 0.8];
      moveAxisShift = [0 0 0 0];
      showTitles = true;
      choiceThreshold = 0.5;
      modelTitleFontSize = 14;
      modelTitleFontName = 'Helvetica';
      modelTitleFontWeight = 'normal';
      titleYNorm = 1.06;
      showXLabel = true;
      showYLabel = true;
      xLabelStr = 'Problem';
      yLabelStr = 'Count LL';
      supAxesPosTemplate = [];
      supAxesBuf = 0.05;
      supYLabelCloser = 0.025; % shift y superlabel right (toward panels)
      supXLabelCloser = 0.025; % shift x superlabel up (toward panels)
      showXTickLabels = true;
      showYTickLabels = true;
      tickLabelsOuterOnly = true;
      sharedYLim = true;
      showGrid = false;
      showSgtitle = false;
      rngSeed = 42;
      regenerateAllModelPostPredSummaries = false;

      summaryPath = fullfile(storageDir, ...
        sprintf('allModelPosteriorPredictive_summary_%s_%s.mat', dataName, engine));

      [pmfPadByModel, matchPropByModel, rebuilt] = loadOrBuildAllModelPostPredSummary( ...
        summaryPath, regenerateAllModelPostPredSummaries, ...
        storageDir, dataName, engine, ...
        cognitiveStems, contaminantStems, cognitiveJagsNames, modelLong, ...
        data, d, nP, nT, nPairs, nModels, nTpMat, nPredSamples, rngSeed, ...
        choiceThreshold);
      if rebuilt
        fprintf('Wrote all-model post-pred summary %s\n', summaryPath);
      else
        fprintf('Loaded all-model post-pred summary %s\n', summaryPath);
      end

      for mi = 1:nModels
        pmfPadByParticipant = pmfPadByModel{mi};
        matchProp = matchPropByModel(mi, :);

        F = figure; clf;
        setFigure(F, figPos, '');
        set(F, 'Color', 'w', 'renderer', 'painters');

        contentPos = nan(nRows * nCols, 4);
        for pp = 1:(nRows * nCols)
          ax = subplot(nRows, nCols, pp);
          if pp > nP
            axis(ax, 'off');
            continue;
          end

          if isempty(pmfPadByParticipant) || isempty(pmfPadByParticipant{pp})
            text(ax, 0.5, 0.5, 'No fit', 'Units', 'normalized', ...
              'HorizontalAlignment', 'center', 'FontSize', fontSize);
            if showTitles
              text(ax, 0, titleYNorm, participantLabels{pp}, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'bottom', ...
                'FontName', titleFontName, ...
                'FontSize', titleFontSize, ...
                'FontWeight', titleFontWeight, ...
                'Interpreter', 'none', ...
                'Clipping', 'off');
            end
            continue;
          end

          pmfPad = pmfPadByParticipant{pp};

          if sharedYLim
            yHi = maxKGlobal;
          else
            yHi = maxKByPar(pp);
          end

          if isempty(pmfPad) || all(isnan(pmfPad(:)), 'all')
            text(ax, 0.5, 0.5, 'No \theta', 'Units', 'normalized', ...
              'HorizontalAlignment', 'center', 'FontSize', fontSize);
          else
            maxProb = max(pmfPad(:), [], 'omitnan');
            if ~(maxProb > 0)
              text(ax, 0.5, 0.5, 'No mass', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', fontSize);
            else
              scale = 0.92 / sqrt(maxProb);
              hold(ax, 'on');

              obsRow = obsCntMat(pp, :);
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
                dispIdx = displayNumber(pj);
                side = min(scale * sqrt(prob), 0.92);
                ko = obsRow(pj);
                isObs = isfinite(ko) && round(ko) == k;
                if isObs
                  continue;
                end
                rectangle(ax, 'Position', [dispIdx - side / 2, k - side / 2, side, side], ...
                  'FaceColor', predFaceClr, 'EdgeColor', predEdgeClr, 'LineWidth', tileLineWidth);
              end

              for dispIdx = 1:nPairs
                pj = problemOrder(dispIdx);
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
                face = problemClr{pairTypeByProblem(pj)};
                edge = max(0, face * (1 - obsEdgeDarken));
                rectangle(ax, 'Position', [dispIdx - side / 2, ko - side / 2, side, side], ...
                  'FaceColor', face, 'EdgeColor', edge, 'LineWidth', obsLineWidth);
              end
              hold(ax, 'off');
            end
          end

          xlim(ax, [1 nPairs]);
          ylim(ax, [0 yHi]);
          set(ax, ...
            'TickDir', 'out', ...
            'Box', 'off', ...
            'FontSize', fontSize, ...
            'XTick', problemXtickLabels, ...
            'xticklabelrot', 0, ...
            'YTick', [0 yHi], ...
            'TickLength', [tickLength 0], ...
            'clipping', 'off');
          onBottom = pp > (nRows - 1) * nCols;
          onLeft = mod(pp, nCols) == 1;
          showXHere = showXTickLabels && (~tickLabelsOuterOnly || onBottom);
          showYHere = showYTickLabels && (~tickLabelsOuterOnly || onLeft);
          if ~showXHere
            set(ax, 'XTickLabel', []);
          end
          if ~showYHere
            set(ax, 'YTickLabel', []);
          end
          if showGrid
            grid(ax, 'on');
          end
          moveAxis(gca, moveAxisScale, moveAxisShift);
          contentPos(pp, :) = get(ax, 'Position');
          [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
          set([axX axY], 'Tag', 'RaxesCopy');
          title(axX, '');
          title(axY, '');
          title(ax, '');
          if showTitles
            text(ax, 0, titleYNorm, participantLabels{pp}, ...
              'Units', 'normalized', ...
              'HorizontalAlignment', 'left', ...
              'VerticalAlignment', 'bottom', ...
              'FontName', titleFontName, ...
              'FontSize', titleFontSize, ...
              'FontWeight', titleFontWeight, ...
              'Interpreter', 'none', ...
              'Clipping', 'off');
            if isfinite(matchProp(pp))
              matchStr = sprintf('%d%%', round(100 * matchProp(pp)));
              text(ax, 1, titleYNorm, matchStr, ...
                'Units', 'normalized', ...
                'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'bottom', ...
                'FontName', modelTitleFontName, ...
                'FontSize', modelTitleFontSize, ...
                'FontWeight', modelTitleFontWeight, ...
                'Interpreter', 'none', ...
                'Clipping', 'off');
            end
          end
        end

        supAxesPos = supAxesPosTemplate;
        if isempty(supAxesPos)
          keep = all(isfinite(contentPos), 2);
          pos = contentPos(keep, :);
          leftMin = min(pos(:, 1));
          bottomMin = min(pos(:, 2));
          leftMax = max(pos(:, 1) + pos(:, 3));
          bottomMax = max(pos(:, 2) + pos(:, 4));
          supAxesPos = [ ...
            leftMin - supAxesBuf - raxesYShift, ...
            bottomMin - supAxesBuf - raxesXShift, ...
            (leftMax - leftMin) + 2 * supAxesBuf + raxesYShift, ...
            (bottomMax - bottomMin) + 2 * supAxesBuf + raxesXShift];
        end

        if showXLabel
          [~, Hx] = suplabel(xLabelStr, 'x', supAxesPos);
          hxPos = get(Hx, 'Position');
          hxPos(2) = hxPos(2) + supXLabelCloser;
          set(Hx, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
            'Position', hxPos);
        end
        if showYLabel
          [~, Hy] = suplabel(yLabelStr, 'y', supAxesPos);
          hyPos = get(Hy, 'Position');
          hyPos(1) = hyPos(1) + supYLabelCloser;
          set(Hy, 'FontSize', labelFontSize, 'VerticalAlignment', 'middle', ...
            'Position', hyPos);
        end

        if showSgtitle
          sgtitle(sprintf('%s posterior predictive counts', modelLong{mi}), ...
            'FontWeight', 'normal', 'FontSize', titleFontSize + 2);
        end

        if printFigures
          pngPath = fullfile(figuresDir, ...
            sprintf('%s_%s.png', analysisName, modelShort{mi}));
          epsPath = fullfile(figuresDir, ...
            sprintf('%s_%s.eps', analysisName, modelShort{mi}));
          print(pngPath, '-dpng', '-r300');
          print(epsPath, '-depsc');
          fprintf('Saved %s and %s\n', pngPath, epsPath);
        end
      end
      saveFigure = false; % already saved per-model above (or interactive)

    %% ================================================================
    case 'modelAgreementTables'
      % APA LaTeX tables (manuscripts/): forced-choice agreement and mean
      % P(observed decision) for every model x participant. Bold = MAP model.
      % Forced choice: mean P(LL) >= .5 => LL, else SS.

      choiceThreshold = 0.5; % >= threshold => LL (so Gu at .5 is determinate)
      manuscriptsDir = fullfile(resultsDir, '..', 'manuscripts');

      mapSummaryPath = fullfile(storageDir, ...
        sprintf('mapModelForcedChoiceMatch_summary_origMAP_original_%s_%s.mat', ...
          dataName, engine));
      if ~isfile(mapSummaryPath)
        mapSummaryPath = fullfile(storageDir, ...
          sprintf('mapModelForcedChoiceMatch_summary_%s_%s_%s.mat', ...
            mixtureName, dataName, engine));
      end
      if ~isfile(mapSummaryPath)
        error('Missing MAP summary (needed to bold MAP cells): %s', mapSummaryPath);
      end
      mapS = load(mapSummaryPath, 'mapModel');
      mapModelRow = mapS.mapModel(:)';
      if numel(mapModelRow) ~= nP
        error('mapModel length %d does not match nP=%d', numel(mapModelRow), nP);
      end
      fprintf('Loaded MAP models from %s\n', mapSummaryPath);

      matchProp = nan(nModels, nP);
      meanObsProb = nan(nModels, nP);
      for mi = 1:nModels
        if mi <= 8
          stem = cognitiveStems{mi};
        else
          stem = contaminantStems{mi - 8};
        end
        fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', stem, dataName, engine));
        if ~isfile(fpath)
          warning('Missing %s', fpath);
          continue;
        end
        fprintf('Loading %s\n', fpath);
        L = load(fpath, 'chains');
        chains = L.chains;
        for pp = 1:nP
          meanTh = mapModelMeanTheta( ...
            chains, data, cognitiveJagsNames, mi, pp, nT);
          if isempty(meanTh) || numel(meanTh) ~= nT
            continue;
          end
          obs = d.LL(pp, :);
          validObs = isfinite(obs) & isfinite(meanTh);

          forcedLL = nan(size(meanTh));
          forcedLL(meanTh >= choiceThreshold) = 1;
          forcedLL(meanTh < choiceThreshold) = 0;
          validForce = isfinite(forcedLL) & validObs;
          if any(validForce)
            matchProp(mi, pp) = mean(forcedLL(validForce) == obs(validForce));
          end

          pObs = nan(size(meanTh));
          pObs(obs == 1) = meanTh(obs == 1);
          pObs(obs == 0) = 1 - meanTh(obs == 0);
          if any(validObs)
            meanObsProb(mi, pp) = mean(pObs(validObs));
          end
        end
        clear chains L;
      end

      abbrevNote = [ ...
        'Model abbreviations: Ex~=~exponential, Hc~=~hyperbolic, Hd~=~hyperboloid, ', ...
        'PD~=~proportional differences, DD~=~direct differences, Tr~=~tradeoff, ', ...
        'UT~=~unified tradeoff, IT~=~intertemporal choice heuristic, ', ...
        'Gu~=~guess, LL~=~larger-later contaminant, SS~=~smaller-sooner contaminant. ', ...
        '\textbf{Bold} entries mark the model with the highest posterior probability ', ...
        'for that participant in the latent-mixture analysis.'];

      writeModelParticipantTable( ...
        fullfile(manuscriptsDir, 'forcedChoiceAgreementTable.tex'), ...
        matchProp, mapModelRow, modelShort, ...
        ['Forced-choice posterior predictive agreement (percentage of trials) ', ...
         'for each model and participant. For each trial, the model''s posterior ', ...
         'mean $P(\mathrm{LL})$ was mapped to a forced larger-later prediction if ', ...
         'it was at least $.5$ and to a smaller-sooner prediction otherwise; ', ...
         'agreement is the proportion of trials matching the participant''s ', ...
         'observed choice. ', abbrevNote], ...
        'tab:forcedChoiceAgreement', true);

      writeModelParticipantTable( ...
        fullfile(manuscriptsDir, 'meanObservedDecisionProbTable.tex'), ...
        meanObsProb, mapModelRow, modelShort, ...
        ['Mean posterior predictive probability of the observed decision ', ...
         '(percentage) for each model and participant. For each trial, the ', ...
         'probability assigned to the participant''s chosen alternative was taken ', ...
         'from the model''s posterior mean $P(\mathrm{LL})$, and these trial-level ', ...
         'probabilities were averaged. ', abbrevNote], ...
        'tab:meanObservedDecisionProb', true);

      writeModelParticipantCsv( ...
        fullfile(figuresDir, 'forcedChoiceAgreement_allModels.csv'), ...
        matchProp, modelShort);
      writeModelParticipantCsv( ...
        fullfile(figuresDir, 'meanObservedDecisionProb_allModels.csv'), ...
        meanObsProb, modelShort);

      saveFigure = false; % tables only

    %% ================================================================
    case 'mapModelForcedChoiceMatch'
      % Forced-choice agreement for each participant's original-MAP model,
      % comparing that same model's fits under original / half / double mu
      % priors (ignores MAP changes under half/double). Cognitive half/double
      % fits: runHierarchicalExecutionPriorRobustness.
      % Console report + two-panel scatter (original vs half / original vs double).

      % graphics / analysis constants
      choiceThreshold = 0.5; % mean theta > threshold => LL; < => SS; == excluded
      % Cutdown summary of match rates (full chain mats kept).
      regenerateMapForcedChoiceSummaries = false;

      mixPathOrig = fullfile(storageDir, sprintf('%s_%s_%s.mat', mixtureName, dataName, engine));
      if ~isfile(mixPathOrig)
        error('Missing mixture chains: %s', mixPathOrig);
      end

      % Stem lists for hierarchical fits under each prior setting
      cognitiveStemsHalf = cellfun(@(s) cognitiveMuPrecStem(s, 'Half'), ...
        cognitiveStems, 'UniformOutput', false);
      cognitiveStemsDouble = cellfun(@(s) cognitiveMuPrecStem(s, 'Double'), ...
        cognitiveStems, 'UniformOutput', false);

      priorPack = { ...
        'Original', 'origMAP_original', cognitiveStems; ...
        'Half',     'origMAP_half',     cognitiveStemsHalf; ...
        'Double',   'origMAP_double',   cognitiveStemsDouble};
      nPrior = size(priorPack, 1);
      matchByMix = cell(nPrior, 1);
      mapModelByMix = cell(nPrior, 1);
      mapProbByMix = cell(nPrior, 1);
      overallByMix = nan(1, nPrior);

      for kPrior = 1:nPrior
        priorLabel = priorPack{kPrior, 1};
        fitTag = priorPack{kPrior, 2};
        cogStems = priorPack{kPrior, 3};
        summaryPath = fullfile(storageDir, ...
          sprintf('mapModelForcedChoiceMatch_summary_%s_%s_%s.mat', ...
            fitTag, dataName, engine));

        [mapModelK, mapProbK, matchPropK, ~, rebuilt] = ...
          loadOrBuildMapForcedChoiceMatchSummary( ...
            summaryPath, regenerateMapForcedChoiceSummaries, ...
            mixPathOrig, storageDir, dataName, engine, ...
            fitTag, cogStems, contaminantStems, cognitiveJagsNames, modelLong, ...
            data, d, nP, nT, nModels, choiceThreshold);
        if rebuilt
          fprintf('Wrote forced-choice match summary (%s) %s\n', priorLabel, summaryPath);
        else
          fprintf('Loaded forced-choice match summary (%s) %s\n', priorLabel, summaryPath);
        end

        fprintf('\nMAP forced-choice match — %s priors (theta > %.2f => LL, < => SS):\n', ...
          priorLabel, choiceThreshold);
        fprintf('  (MAP identity from original mixture; params from %s hierarchical fits)\n', ...
          lower(priorLabel));
        for pp = 1:nP
          mi = mapModelK(pp);
          if isfinite(mi) && mi >= 1 && mi <= nModels && isfinite(matchPropK(pp))
            fprintf('  %s  %s[%.2f]  match=%.3f\n', ...
              participantLabels{pp}, modelShort{mi}, mapProbK(pp), matchPropK(pp));
          else
            fprintf('  %s  (no MAP fit)\n', participantLabels{pp});
          end
        end
        overallK = mean(matchPropK(isfinite(matchPropK)));
        nReady = sum(isfinite(matchPropK));
        fprintf('  Overall mean match = %.3f (%d / %d participants with fits ready)\n', ...
          overallK, nReady, nP);

        matchByMix{kPrior} = matchPropK;
        mapModelByMix{kPrior} = mapModelK;
        mapProbByMix{kPrior} = mapProbK;
        overallByMix(kPrior) = overallK;
      end

      fprintf('\nForced-choice agreement summary (same original-MAP model across priors):\n');
      for kPrior = 1:nPrior
        fprintf('  %-8s overall mean match = %.3f\n', ...
          [priorPack{kPrior, 1} ':'], overallByMix(kPrior));
      end
      fprintf('\n');

      % Two-panel scatter: original vs half / original vs double (% agreement)
      fontSize = 14;
      labelFontSize = 16;
      titleFontSize = 14;
      figPos = [0.15 0.25 0.7 0.45];
      markerSize = 36;
      markerFace = [0.45 0.55 0.75];
      markerEdge = 0.25 * [1 1 1];
      identityClr = 0.55 * [1 1 1];
      identityLineWidth = 1.0;
      tickLength = 0.02;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 0.9 0.9];
      moveAxisShift = [0 0.025 0 0];
      xLabelStr = 'Original prior';
      yLabelHalf = 'Half prior';
      yLabelDouble = 'Double prior';
      showPanelTitles = true;
      panelTitleHalf = 'Half';
      panelTitleDouble = 'Double';

      pctOrig = 100 * matchByMix{1};
      pctHalf = 100 * matchByMix{2};
      pctDouble = 100 * matchByMix{3};
      altPack = {pctHalf, pctDouble};
      yLabs = {yLabelHalf, yLabelDouble};
      panTitles = {panelTitleHalf, panelTitleDouble};

      F = figure; clf;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      for pan = 1:2
        ax = subplot(1, 2, pan);
        hold(ax, 'on');
        pctAlt = altPack{pan};

        plot(ax, [0 100], [0 100], '-', 'Color', identityClr, 'LineWidth', identityLineWidth);

        for pp = 1:nP
          xo = pctOrig(pp);
          ya = pctAlt(pp);
          if ~(isfinite(xo) && isfinite(ya))
            continue;
          end
          scatter(ax, xo, ya, markerSize, ...
            'MarkerFaceColor', markerFace, ...
            'MarkerEdgeColor', markerEdge, ...
            'LineWidth', 0.6);
        end

        xlim(ax, [0 100]);
        ylim(ax, [0 100]);
        axis(ax, 'square');
        set(ax, ...
          'TickDir', 'out', ...
          'Box', 'off', ...
          'FontSize', fontSize, ...
          'XTick', 0:25:100, ...
          'YTick', 0:25:100, ...
          'TickLength', [tickLength 0], ...
          'clipping', 'off');
        xlabel(ax, xLabelStr, 'FontSize', labelFontSize);
        ylabel(ax, yLabs{pan}, 'FontSize', labelFontSize);
        if showPanelTitles
          title(ax, panTitles{pan}, 'FontSize', titleFontSize, 'FontWeight', 'normal');
        end
        moveAxis(gca, moveAxisScale, moveAxisShift);
        [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
        set([axX axY], 'Tag', 'RaxesCopy');
        title(axX, '');
        title(axY, '');
        if showPanelTitles
          title(ax, panTitles{pan}, 'FontSize', titleFontSize, 'FontWeight', 'normal');
        end
      end

    %% ================================================================
    case 'mapModelParameterRobustness'
      % For each parameter of each model that is MAP for someone under the
      % original prior: two-panel scatter of participant posterior means
      % (original vs half / original vs double), only those MAP participants.
      % Always uses that original-MAP model for all three prior settings
      % (ignores MAP changes under half/double). Cognitive half/double fits
      % come from runHierarchicalExecutionPriorRobustness.

      fontSize = 9;
      labelFontSize = 9;
      titleFontSize = 9;
      figPos = [0.02 0.02 0.96 0.96];
      markerSize = 28;
      nHalfCols = 3;
      nDoubleCols = 3;
      markerFace = [0.45 0.55 0.75];
      markerEdge = 0.25 * [1 1 1];
      identityClr = 0.55 * [1 1 1];
      identityLineWidth = 0.8;
      tickLength = 0.02;
      raxesXShift = 0.008;
      raxesYShift = 0.008;
      moveAxisScale = [1 1 1 1];
      moveAxisShift = [0 0 0 0];
      rhsLabelXNorm = 1.12; % normalized x for right-hand panel title
      xLabelStr = 'Original prior';
      yLabelHalf = 'Half prior';
      yLabelDouble = 'Double prior';
      errBarClr = 0.65 * [1 1 1];
      errBarLineWidth = 0.7;
      axisPadFrac = 0.08;
      iqrProbs = [0.25 0.75];

      % {label, individualFitField, mixtureField} per cognitive/contaminant model
      modelParamSpecs = { ...
        {{'kappa', 'kappa', 'kappaEX'}, {'w', 'w', 'w'}}; ... % Ex
        {{'kappa', 'kappa', 'kappaHC'}, {'w', 'w', 'w'}}; ... % Hc
        {{'kappa', 'kappa', 'kappaHY'}, {'tau', 'tau', 'tauHY'}, {'w', 'w', 'w'}}; ... % Hd
        {{'delta', 'delta', 'deltaPD'}, {'w', 'w', 'w'}}; ... % PD
        {{'delta', 'delta', 'deltaDD'}, {'omega', 'omega', 'omegaDD'}, {'w', 'w', 'w'}}; ... % DD
        {{'gamma', 'gamma', 'gammaTR'}, {'tau', 'tau', 'tauTR'}, ...
          {'kappa', 'kappa', 'kappaTR'}, {'vartheta', 'vartheta', 'varthetaTR'}, {'w', 'w', 'w'}}; ... % Tr
        {{'gamma', 'gamma', 'gammaUT'}, {'tau', 'tau', 'tauUT'}, ...
          {'kappa', 'kappa', 'kappaUT'}, {'vartheta', 'vartheta', 'varthetaUT'}, ...
          {'eta', 'eta', 'etaUT'}, {'w', 'w', 'w'}}; ... % UT
        {{'beta0', 'beta0', 'beta0IT'}, {'betaRA', 'betaRA', 'betaRAIT'}, ...
          {'betaRR', 'betaRR', 'betaRRIT'}, {'betaTA', 'betaTA', 'betaTAIT'}, ...
          {'betaTR', 'betaTR', 'betaTRIT'}, {'w', 'w', 'w'}}; ... % IT
        {}; ... % Gu
        {{'alpha', 'alpha', 'alphaLL'}}; ... % LL
        {{'alpha', 'alpha', 'alphaSS'}} ... % SS
        };

      % MAP model under original mixture priors only
      mixPath = fullfile(storageDir, sprintf('%s_%s_%s.mat', mixtureName, dataName, engine));
      if ~isfile(mixPath)
        error('Missing mixture chains: %s', mixPath);
      end
      fprintf('Loading mixture %s\n', mixPath);
      mixS = load(mixPath, 'chains');
      P = posteriorModelProbsFromZ(mixS.chains, nP, nModels);
      [~, mapOrig] = max(P, [], 1);
      clear mixS P;

      uniqueMaps = unique(mapOrig(isfinite(mapOrig)));
      uniqueMaps = uniqueMaps(uniqueMaps >= 1 & uniqueMaps <= nModels);
      if isempty(uniqueMaps)
        error('No MAP models found under original mixture');
      end
      fprintf('Original MAP models: %s\n', strjoin(lower(modelShort(uniqueMaps)), ', '));

      % Load original / half / double hierarchical fits for each needed MAP model.
      % Half/double cognitive fits: runHierarchicalExecutionPriorRobustness.
      % Skip models whose half/double fits are not on disk yet (incremental runs).
      chainsOrig = cell(nModels, 1);
      chainsHalf = cell(nModels, 1);
      chainsDouble = cell(nModels, 1);
      for mi = uniqueMaps(:)'
        specs = modelParamSpecs{mi};
        if isempty(specs)
          continue;
        end
        if mi > 8
          fprintf('Skipping contaminant MAP model %s (no mu_* prior-robustness fits)\n', ...
            modelShort{mi});
          continue;
        end
        stemOrig = cognitiveStems{mi};
        stemHalf = cognitiveMuPrecStem(stemOrig, 'Half');
        stemDouble = cognitiveMuPrecStem(stemOrig, 'Double');
        hasHalf = hierarchicalStemAvailable(storageDir, stemHalf, dataName, engine);
        hasDouble = hierarchicalStemAvailable(storageDir, stemDouble, dataName, engine);
        if ~(hasHalf || hasDouble)
          fprintf('Skipping %s — no half/double fit yet\n', modelShort{mi});
          continue;
        end
        if ~hierarchicalStemAvailable(storageDir, stemOrig, dataName, engine)
          fprintf('Skipping %s — missing original fit\n', modelShort{mi});
          continue;
        end
        chainsOrig{mi} = loadHierarchicalChains(storageDir, stemOrig, dataName, engine);
        if hasHalf
          chainsHalf{mi} = loadHierarchicalChains(storageDir, stemHalf, dataName, engine);
        else
          fprintf('  %s: half fit not ready — half panel points will be omitted\n', ...
            modelShort{mi});
        end
        if hasDouble
          chainsDouble{mi} = loadHierarchicalChains(storageDir, stemDouble, dataName, engine);
        else
          fprintf('  %s: double fit not ready — double panel points will be omitted\n', ...
            modelShort{mi});
        end
      end

      % Build list of (model, param) panels to plot
      panelJobs = {};
      for mi = uniqueMaps(:)'
        specs = modelParamSpecs{mi};
        if isempty(specs) || mi > 8
          if isempty(specs)
            fprintf('Skipping model %s (no plotted parameters)\n', modelShort{mi});
          end
          continue;
        end
        if isempty(chainsOrig{mi})
          continue;
        end
        hasHalf = ~isempty(chainsHalf{mi});
        hasDouble = ~isempty(chainsDouble{mi});
        if ~(hasHalf || hasDouble)
          continue;
        end
        pps = find(mapOrig == mi);
        if isempty(pps)
          continue;
        end
        for s = 1:numel(specs)
          label = specs{s}{1};
          indivField = specs{s}{2};
          muOrig = nan(1, nP);
          loOrig = nan(1, nP);
          hiOrig = nan(1, nP);
          muHalf = nan(1, nP);
          loHalf = nan(1, nP);
          hiHalf = nan(1, nP);
          muDouble = nan(1, nP);
          loDouble = nan(1, nP);
          hiDouble = nan(1, nP);
          for ii = 1:numel(pps)
            pp = pps(ii);
            [muOrig(pp), loOrig(pp), hiOrig(pp)] = indexedParamSummary( ...
              chainsOrig{mi}, indivField, pp, iqrProbs);
            if hasHalf
              [muHalf(pp), loHalf(pp), hiHalf(pp)] = indexedParamSummary( ...
                chainsHalf{mi}, indivField, pp, iqrProbs);
            end
            if hasDouble
              [muDouble(pp), loDouble(pp), hiDouble(pp)] = indexedParamSummary( ...
                chainsDouble{mi}, indivField, pp, iqrProbs);
            end
          end
          panelJobs{end + 1} = struct( ... %#ok<AGROW>
            'mi', mi, ...
            'label', label, ...
            'pps', pps, ...
            'muOrig', muOrig, 'loOrig', loOrig, 'hiOrig', hiOrig, ...
            'muHalf', muHalf, 'loHalf', loHalf, 'hiHalf', hiHalf, ...
            'muDouble', muDouble, 'loDouble', loDouble, 'hiDouble', hiDouble);
        end
      end

      nJobs = numel(panelJobs);
      if nJobs == 0
        fprintf(['No parameter panels to plot yet (waiting on ' ...
          'runHierarchicalExecutionPriorRobustness fits).\n']);
        continue;
      end
      nFigCols = nHalfCols + nDoubleCols;
      nFigRows = ceil(nJobs / nHalfCols);
      fprintf('Drawing %d parameters on a %d x %d grid (3 half + 3 double cols)\n', ...
        nJobs, nFigRows, nFigCols);

      F = figure; clf;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      for j = 1:nJobs
        job = panelJobs{j};
        row = ceil(j / nHalfCols);
        colWithin = mod(j - 1, nHalfCols) + 1;
        altMu = {job.muHalf, job.muDouble};
        altLo = {job.loHalf, job.loDouble};
        altHi = {job.hiHalf, job.hiDouble};
        yLabs = {yLabelHalf, yLabelDouble};
        subplotIdx = [(row - 1) * nFigCols + colWithin, ...
                      (row - 1) * nFigCols + nHalfCols + colWithin];
        for pan = 1:2
          ax = subplot(nFigRows, nFigCols, subplotIdx(pan));
          hold(ax, 'on');
          muAlt = altMu{pan};
          loAlt = altLo{pan};
          hiAlt = altHi{pan};
          xs = [];
          ys = [];
          for ii = 1:numel(job.pps)
            pp = job.pps(ii);
            xo = job.muOrig(pp);
            ya = muAlt(pp);
            if ~(isfinite(xo) && isfinite(ya))
              continue;
            end
            % IQR error bars from marginal posteriors (x = original, y = alt)
            xLo = job.loOrig(pp);
            xHi = job.hiOrig(pp);
            yLo = loAlt(pp);
            yHi = hiAlt(pp);
            if isfinite(xLo) && isfinite(xHi)
              plot(ax, [xLo xHi], [ya ya], '-', ...
                'Color', errBarClr, 'LineWidth', errBarLineWidth);
            end
            if isfinite(yLo) && isfinite(yHi)
              plot(ax, [xo xo], [yLo yHi], '-', ...
                'Color', errBarClr, 'LineWidth', errBarLineWidth);
            end
            scatter(ax, xo, ya, markerSize, ...
              'MarkerFaceColor', markerFace, ...
              'MarkerEdgeColor', markerEdge, ...
              'LineWidth', 0.5);
            xs(end + 1) = xo; %#ok<AGROW>
            ys(end + 1) = ya; %#ok<AGROW>
            if isfinite(xLo), xs(end + 1) = xLo; end %#ok<AGROW>
            if isfinite(xHi), xs(end + 1) = xHi; end %#ok<AGROW>
            if isfinite(yLo), ys(end + 1) = yLo; end %#ok<AGROW>
            if isfinite(yHi), ys(end + 1) = yHi; end %#ok<AGROW>
          end

          if isempty(xs)
            text(ax, 0.5, 0.5, 'No paired means', 'Units', 'normalized', ...
              'HorizontalAlignment', 'center', 'FontSize', fontSize);
            xlim(ax, [0 1]);
            ylim(ax, [0 1]);
          else
            lo = min([xs, ys]);
            hi = max([xs, ys]);
            if ~(isfinite(lo) && isfinite(hi)) || lo == hi
              pad = max(abs(lo), 1);
              lo = lo - pad;
              hi = hi + pad;
            else
              pad = axisPadFrac * (hi - lo);
              lo = lo - pad;
              hi = hi + pad;
            end
            plot(ax, [lo hi], [lo hi], '-', 'Color', identityClr, ...
              'LineWidth', identityLineWidth);
            xlim(ax, [lo hi]);
            ylim(ax, [lo hi]);
          end
          axis(ax, 'square');
          set(ax, ...
            'TickDir', 'out', ...
            'Box', 'off', ...
            'FontSize', fontSize, ...
            'TickLength', [tickLength 0], ...
            'clipping', 'off');
          title(ax, '');
          if row == nFigRows
            xlabel(ax, xLabelStr, 'FontSize', labelFontSize);
          else
            xlabel(ax, '');
          end
          ylabel(ax, yLabs{pan}, 'FontSize', labelFontSize);
          moveAxis(gca, moveAxisScale, moveAxisShift);
          [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
          set([axX axY], 'Tag', 'RaxesCopy');
          title(axX, '');
          title(axY, '');
          title(ax, '');
          % Panel title as right-hand y-axis label
          text(ax, rhsLabelXNorm, 0.5, ...
            sprintf('%s %s (n=%d)', lower(modelShort{job.mi}), job.label, numel(job.pps)), ...
            'Units', 'normalized', ...
            'Rotation', -90, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', titleFontSize, ...
            'FontWeight', 'normal', ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
        end
      end

    %% ================================================================
    case 'parameterInferences'
      % Per-participant MAP-model parameter means (+ 95% CI), adapted from
      % drawFigures_4 parameterInferences, with entropification w (not epsilon).

      fontSize = 10;
      tickLabelFontSize = fontSize + 1; % LaTeX ticks need a bit more size
      xTickLabelRotation = 90;
      titleFontSize = 14;
      titleFontName = 'Helvetica';
      titleFontWeight = 'normal';
      modelTitleFontSize = 14;
      modelTitleFontName = 'Helvetica';
      modelTitleFontWeight = 'normal';
      titleYNorm = 1.06;
      figPos = [0.2 0.2 0.45 0.6];
      CIbounds = [2.5 97.5];
      markerFace = pantone.ClassicBlue;
      markerSize = 5;
      ciLineWidth = 1.2;
      tickLength = 0.02;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 1 0.7];
      moveAxisShift = [0 0.015 0 0];
      yLimFixed = [-5 10]; % match drawFigures_4; set [] for per-panel auto
      meanLabelFontSize = 7;
      meanLabelYOffset = 0.45; % data units above the mean marker
      nRows = 5;
      nCols = 5;

      % Fields in hierarchical / contaminant chain mats (epsilon -> w)
      parameterNames = { ...
        {'kappa', 'w'}; ...                                  % Ex
        {'kappa', 'w'}; ...                                  % Hc
        {'kappa', 'tau', 'w'}; ...                           % Hd
        {'delta', 'w'}; ...                                  % PD
        {'delta', 'omega', 'w'}; ...                         % DD
        {'gamma', 'kappa', 'tau', 'vartheta', 'w'}; ...       % Tr
        {'gamma', 'kappa', 'tau', 'vartheta', 'eta', 'w'}; ...% UT
        {'betaRA', 'betaRR', 'betaTA', 'betaTR', 'beta0', 'w'}; ... % IT
        {}; ...                                              % Gu
        {'alpha'}; ...                                       % LL
        {'alpha'} ...                                        % SS
        };

      modelStems = [cognitiveStems, contaminantStems];
      mn = cell(nModels, 1);
      ci = cell(nModels, 1);
      for mi = 1:nModels
        nPar = numel(parameterNames{mi});
        mn{mi} = nan(nP, max(nPar, 1));
        ci{mi} = nan(nP, max(nPar, 1), 2);
        if nPar == 0
          continue;
        end
        fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', modelStems{mi}, dataName, engine));
        if ~isfile(fpath)
          error('Missing model chains for %s: %s', modelShort{mi}, fpath);
        end
        fprintf('Loading %s\n', fpath);
        S = load(fpath, 'chains');
        for pp = 1:nP
          for k = 1:nPar
            [mu, qLo, qHi] = indexedParamSummary( ...
              S.chains, parameterNames{mi}{k}, pp, CIbounds / 100);
            mn{mi}(pp, k) = mu;
            ci{mi}(pp, k, 1) = qLo;
            ci{mi}(pp, k, 2) = qHi;
          end
        end
        clear S;
      end

      mixPath = fullfile(storageDir, sprintf('%s_%s_%s.mat', mixtureName, dataName, engine));
      if ~isfile(mixPath)
        error('Missing mixture chains: %s', mixPath);
      end
      fprintf('Loading MAP labels from %s\n', mixPath);
      mixS = load(mixPath, 'chains');
      P = posteriorModelProbsFromZ(mixS.chains, nP, nModels);
      [~, mapModel] = max(P, [], 1);
      clear mixS P;

      F = figure; clf;
      setFigure(F, figPos, '');
      set(F, 'Color', 'w', 'renderer', 'painters');

      for pp = 1:(nRows * nCols)
        ax = subplot(nRows, nCols, pp);
        if pp > nP
          axis(ax, 'off');
          continue;
        end
        mi = mapModel(pp);
        pNames = parameterNames{mi};
        nPar = numel(pNames);
        hold(ax, 'on');
        if nPar == 0
          text(ax, 0.5, 0.5, lower(modelShort{mi}), ...
            'Units', 'normalized', 'HorizontalAlignment', 'center', ...
            'FontSize', titleFontSize);
          set(ax, 'XTick', [], 'YTick', [], 'Box', 'off');
          hold(ax, 'off');
          text(ax, 0, titleYNorm, participantLabels{pp}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', titleFontName, ...
            'FontSize', titleFontSize, ...
            'FontWeight', titleFontWeight, ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
          text(ax, 1, titleYNorm, lower(modelShort{mi}), ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', modelTitleFontName, ...
            'FontSize', modelTitleFontSize, ...
            'FontWeight', modelTitleFontWeight, ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
          continue;
        end
        for k = 1:nPar
          plot(ax, [k k], squeeze(ci{mi}(pp, k, :))', '-', ...
            'Color', markerFace, 'LineWidth', ciLineWidth);
          plot(ax, k, mn{mi}(pp, k), 'o', ...
            'MarkerFaceColor', markerFace, ...
            'MarkerEdgeColor', 'none', ...
            'MarkerSize', markerSize);
          if isfinite(mn{mi}(pp, k))
            text(ax, k, mn{mi}(pp, k) + meanLabelYOffset, ...
              sprintf('%.2g', mn{mi}(pp, k)), ...
              'HorizontalAlignment', 'center', ...
              'VerticalAlignment', 'bottom', ...
              'FontSize', meanLabelFontSize, ...
              'Color', markerFace, ...
              'Clipping', 'off');
          end
        end
        hold(ax, 'off');
        xlim(ax, [1 nPar] + [-0.5 0.5]);
        if isempty(yLimFixed)
          vals = [mn{mi}(pp, 1:nPar), reshape(ci{mi}(pp, 1:nPar, :), 1, [])];
          vals = vals(isfinite(vals));
          if isempty(vals)
            yl = [-1 1];
          else
            pad = 0.1 * max(range(vals), 1);
            yl = [min(vals) - pad, max(vals) + pad];
          end
          ylim(ax, yl);
        else
          ylim(ax, yLimFixed);
        end
        set(ax, ...
          'XTick', 1:nPar, ...
          'XTickLabel', parameterLatexLabels(pNames), ...
          'XTickLabelRotation', xTickLabelRotation, ...
          'TickLabelInterpreter', 'latex', ...
          'TickDir', 'out', ...
          'TickLength', [tickLength 0], ...
          'Box', 'off', ...
          'FontSize', tickLabelFontSize, ...
          'Clipping', 'off');
        title(ax, '');
        moveAxis(gca, moveAxisScale, moveAxisShift);
        [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
        set([axX axY], 'Tag', 'RaxesCopy');
        set(axX, 'TickLabelInterpreter', 'latex', ...
          'XTickLabelRotation', xTickLabelRotation, ...
          'FontSize', tickLabelFontSize);
        title(axX, '');
        title(axY, '');
        title(ax, '');
        text(ax, 0, titleYNorm, participantLabels{pp}, ...
          'Units', 'normalized', ...
          'HorizontalAlignment', 'left', ...
          'VerticalAlignment', 'bottom', ...
          'FontName', titleFontName, ...
          'FontSize', titleFontSize, ...
          'FontWeight', titleFontWeight, ...
          'Interpreter', 'none', ...
          'Clipping', 'off');
        text(ax, 1, titleYNorm, lower(modelShort{mi}), ...
          'Units', 'normalized', ...
          'HorizontalAlignment', 'right', ...
          'VerticalAlignment', 'bottom', ...
          'FontName', modelTitleFontName, ...
          'FontSize', modelTitleFontSize, ...
          'FontWeight', modelTitleFontWeight, ...
          'Interpreter', 'none', ...
          'Clipping', 'off');
      end

    %% ================================================================
    case 'allParameterInferences'
      % One 5x5 figure per model: that model's parameter means (+ 95% CI)
      % for every participant (parallel to allModelPosteriorPredictive).

      fontSize = 10;
      tickLabelFontSize = fontSize + 1;
      xTickLabelRotation = 90;
      titleFontSize = 14;
      titleFontName = 'Helvetica';
      titleFontWeight = 'normal';
      modelTitleFontSize = 14;
      modelTitleFontName = 'Helvetica';
      modelTitleFontWeight = 'normal';
      titleYNorm = 1.06;
      figPos = [0.2 0.2 0.45 0.6];
      CIbounds = [2.5 97.5];
      markerFace = pantone.ClassicBlue;
      markerSize = 5;
      ciLineWidth = 1.2;
      tickLength = 0.02;
      raxesXShift = 0.01;
      raxesYShift = 0.01;
      moveAxisScale = [1 1 1 0.7];
      moveAxisShift = [0 0.015 0 0];
      yLimFixed = [-5 10];
      meanLabelFontSize = 7;
      meanLabelYOffset = 0.45;
      nRows = 5;
      nCols = 5;

      parameterNames = { ...
        {'kappa', 'w'}; ...
        {'kappa', 'w'}; ...
        {'kappa', 'tau', 'w'}; ...
        {'delta', 'w'}; ...
        {'delta', 'omega', 'w'}; ...
        {'gamma', 'kappa', 'tau', 'vartheta', 'w'}; ...
        {'gamma', 'kappa', 'tau', 'vartheta', 'eta', 'w'}; ...
        {'betaRA', 'betaRR', 'betaTA', 'betaTR', 'beta0', 'w'}; ...
        {}; ...
        {'alpha'}; ...
        {'alpha'} ...
        };

      modelStems = [cognitiveStems, contaminantStems];
      mn = cell(nModels, 1);
      ci = cell(nModels, 1);
      for mi = 1:nModels
        nPar = numel(parameterNames{mi});
        mn{mi} = nan(nP, max(nPar, 1));
        ci{mi} = nan(nP, max(nPar, 1), 2);
        if nPar == 0
          continue;
        end
        fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', modelStems{mi}, dataName, engine));
        if ~isfile(fpath)
          error('Missing model chains for %s: %s', modelShort{mi}, fpath);
        end
        fprintf('Loading %s\n', fpath);
        S = load(fpath, 'chains');
        for pp = 1:nP
          for k = 1:nPar
            [mu, qLo, qHi] = indexedParamSummary( ...
              S.chains, parameterNames{mi}{k}, pp, CIbounds / 100);
            mn{mi}(pp, k) = mu;
            ci{mi}(pp, k, 1) = qLo;
            ci{mi}(pp, k, 2) = qHi;
          end
        end
        clear S;
      end

      for mi = 1:nModels
        pNames = parameterNames{mi};
        nPar = numel(pNames);

        F = figure; clf;
        setFigure(F, figPos, '');
        set(F, 'Color', 'w', 'renderer', 'painters');

        for pp = 1:(nRows * nCols)
          ax = subplot(nRows, nCols, pp);
          if pp > nP
            axis(ax, 'off');
            continue;
          end
          hold(ax, 'on');
          if nPar == 0
            text(ax, 0.5, 0.5, lower(modelShort{mi}), ...
              'Units', 'normalized', 'HorizontalAlignment', 'center', ...
              'FontSize', titleFontSize);
            set(ax, 'XTick', [], 'YTick', [], 'Box', 'off');
            hold(ax, 'off');
            text(ax, 0, titleYNorm, participantLabels{pp}, ...
              'Units', 'normalized', ...
              'HorizontalAlignment', 'left', ...
              'VerticalAlignment', 'bottom', ...
              'FontName', titleFontName, ...
              'FontSize', titleFontSize, ...
              'FontWeight', titleFontWeight, ...
              'Interpreter', 'none', ...
              'Clipping', 'off');
            continue;
          end
          for k = 1:nPar
            plot(ax, [k k], squeeze(ci{mi}(pp, k, :))', '-', ...
              'Color', markerFace, 'LineWidth', ciLineWidth);
            plot(ax, k, mn{mi}(pp, k), 'o', ...
              'MarkerFaceColor', markerFace, ...
              'MarkerEdgeColor', 'none', ...
              'MarkerSize', markerSize);
            if isfinite(mn{mi}(pp, k))
              text(ax, k, mn{mi}(pp, k) + meanLabelYOffset, ...
                sprintf('%.2g', mn{mi}(pp, k)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', meanLabelFontSize, ...
                'Color', markerFace, ...
                'Clipping', 'off');
            end
          end
          hold(ax, 'off');
          xlim(ax, [1 nPar] + [-0.5 0.5]);
          if isempty(yLimFixed)
            vals = [mn{mi}(pp, 1:nPar), reshape(ci{mi}(pp, 1:nPar, :), 1, [])];
            vals = vals(isfinite(vals));
            if isempty(vals)
              yl = [-1 1];
            else
              pad = 0.1 * max(range(vals), 1);
              yl = [min(vals) - pad, max(vals) + pad];
            end
            ylim(ax, yl);
          else
            ylim(ax, yLimFixed);
          end
          set(ax, ...
            'XTick', 1:nPar, ...
            'XTickLabel', parameterLatexLabels(pNames), ...
            'XTickLabelRotation', xTickLabelRotation, ...
            'TickLabelInterpreter', 'latex', ...
            'TickDir', 'out', ...
            'TickLength', [tickLength 0], ...
            'Box', 'off', ...
            'FontSize', tickLabelFontSize, ...
            'Clipping', 'off');
          title(ax, '');
          moveAxis(gca, moveAxisScale, moveAxisShift);
          [axX, axY] = Raxes(ax, raxesXShift, raxesYShift);
          set([axX axY], 'Tag', 'RaxesCopy');
          set(axX, 'TickLabelInterpreter', 'latex', ...
            'XTickLabelRotation', xTickLabelRotation, ...
            'FontSize', tickLabelFontSize);
          title(axX, '');
          title(axY, '');
          title(ax, '');
          text(ax, 0, titleYNorm, participantLabels{pp}, ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'bottom', ...
            'FontName', titleFontName, ...
            'FontSize', titleFontSize, ...
            'FontWeight', titleFontWeight, ...
            'Interpreter', 'none', ...
            'Clipping', 'off');
        end

        if printFigures
          pngPath = fullfile(figuresDir, ...
            sprintf('%s_%s.png', analysisName, modelShort{mi}));
          epsPath = fullfile(figuresDir, ...
            sprintf('%s_%s.eps', analysisName, modelShort{mi}));
          print(pngPath, '-dpng', '-r300');
          print(epsPath, '-depsc');
          fprintf('Saved %s and %s\n', pngPath, epsPath);
        end
      end
      saveFigure = false; % already saved per-model above

  end

  % ---- print ----
  if printFigures && saveFigure && ~isempty(get(groot, 'CurrentFigure'))
    pngPath = fullfile(figuresDir, sprintf('%s.png', analysisName));
    epsPath = fullfile(figuresDir, sprintf('%s.eps', analysisName));
    print(pngPath, '-dpng', '-r300');
    print(epsPath, '-depsc');
    fprintf('Saved %s and %s\n', pngPath, epsPath);
  end

end

if ~(exist('drawFiguresEntropKeepWorkspace', 'var') && ...
    logical(drawFiguresEntropKeepWorkspace))
  rmpath(generalDir);
  rmpath(modelsDir);
end

%% ----- local helpers -----
function labels = parameterLatexLabels(paramNames)
%PARAMETERLATEXLABELS  Map chain field names to LaTeX math tick labels.
labels = cell(size(paramNames));
for k = 1:numel(paramNames)
  switch paramNames{k}
    case 'kappa',    labels{k} = '$\kappa$';
    case 'tau',      labels{k} = '$\tau$';
    case 'delta',    labels{k} = '$\delta$';
    case 'omega',    labels{k} = '$\omega$';
    case 'gamma',    labels{k} = '$\gamma$';
    case 'vartheta', labels{k} = '$\vartheta$';
    case 'eta',      labels{k} = '$\eta$';
    case 'betaRA',   labels{k} = '$\beta_{RA}$';
    case 'betaRR',   labels{k} = '$\beta_{RR}$';
    case 'betaTA',   labels{k} = '$\beta_{TA}$';
    case 'betaTR',   labels{k} = '$\beta_{TR}$';
    case 'beta0',    labels{k} = '$\beta_0$';
    case 'alpha',    labels{k} = '$\alpha$';
    case 'w',        labels{k} = '$w$';
    otherwise,       labels{k} = sprintf('$\\mathrm{%s}$', paramNames{k});
  end
end
end

function [mapModel, mapProb, pmfPadByParticipant, rebuilt] = loadOrBuildMapPostPredSummary( ...
  summaryPath, forceRegen, ...
  mixPath, storageDir, dataName, engine, ...
  mixtureName, cognitiveStems, contaminantStems, cognitiveJagsNames, modelLong, ...
  data, d, nP, nT, nPairs, nModels, nTpMat, nPredSamples, rngSeed)
%LOADORBUILDMAPPOSTPREDSUMMARY  Cache MAP labels + count PMFs for fast replotting.

needBuild = forceRegen || ~isfile(summaryPath);

if ~needBuild
  S = load(summaryPath);
  needBuild = ~(isfield(S, 'mapModel') && isfield(S, 'mapProb') && ...
    isfield(S, 'pmfPadByParticipant') && isfield(S, 'sourcePaths') && ...
    isfield(S, 'sourceMtimes') && isfield(S, 'nPredSamples') && ...
    isfield(S, 'rngSeed') && isfield(S, 'nP') && isfield(S, 'mixtureName'));
  if ~needBuild
    needBuild = ~(S.nPredSamples == nPredSamples && S.rngSeed == rngSeed && ...
      S.nP == nP && strcmp(S.mixtureName, mixtureName) && ...
      ~isempty(S.sourcePaths) && strcmp(S.sourcePaths{1}, mixPath));
  end
  if ~needBuild
    curMt = nan(numel(S.sourcePaths), 1);
    for k = 1:numel(S.sourcePaths)
      curMt(k) = fileMtime(S.sourcePaths{k});
    end
    needBuild = any(~isfinite(curMt)) || ~isequal(curMt(:), S.sourceMtimes(:));
  end
  if ~needBuild
    mapModel = S.mapModel;
    mapProb = S.mapProb;
    pmfPadByParticipant = S.pmfPadByParticipant;
    rebuilt = false;
    return;
  end
end

fprintf('Building MAP post-pred summary from full chains...\n');
rng(rngSeed);

fprintf('Loading mixture %s\n', mixPath);
mixS = load(mixPath, 'chains');
postModelProb = posteriorModelProbsFromZ(mixS.chains, nP, nModels);
[mapProb, mapModel] = max(postModelProb, [], 1);
clear mixS;

uniqueMaps = unique(mapModel(~isnan(mapModel)));
chainsByModel = cell(nModels, 1);
modelPaths = {};
modelMtimes = [];
for mi = uniqueMaps(:)'
  if mi <= 8
    stem = cognitiveStems{mi};
  else
    stem = contaminantStems{mi - 8};
  end
  fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', stem, dataName, engine));
  if ~isfile(fpath)
    warning('Missing fit for MAP model %s (%s)', modelLong{mi}, fpath);
    continue;
  end
  fprintf('Loading %s\n', fpath);
  L = load(fpath, 'chains');
  chainsByModel{mi} = L.chains;
  modelPaths{end + 1} = fpath; %#ok<AGROW>
  modelMtimes(end + 1, 1) = fileMtime(fpath); %#ok<AGROW>
end

pmfPadByParticipant = cell(1, nP);
for pp = 1:nP
  mi = mapModel(pp);
  if ~(isfinite(mi) && mi >= 1 && mi <= nModels) || isempty(chainsByModel{mi})
    pmfPadByParticipant{pp} = [];
    continue;
  end
  if mi <= 8
    Th = hierarchicalExecutionThetaDraws( ...
      chainsByModel{mi}, data, cognitiveJagsNames{mi}, pp);
    pmfPadByParticipant{pp} = postPredictiveCountPmfFromTheta( ...
      Th, d.pair(pp, :), nPairs, nTpMat(pp, :), nPredSamples);
  else
    pmfPadByParticipant{pp} = postPredictiveCountPmfFromMonitoredTheta( ...
      chainsByModel{mi}, pp, d.pair(pp, :), nT, nPairs, nTpMat(pp, :), nPredSamples);
  end
end
clear chainsByModel;

sourcePaths = [{mixPath}; modelPaths(:)];
sourceMtimes = [fileMtime(mixPath); modelMtimes(:)];
save(summaryPath, ...
  'mapModel', 'mapProb', 'pmfPadByParticipant', ...
  'sourcePaths', 'sourceMtimes', ...
  'nPredSamples', 'rngSeed', 'nP', 'mixtureName', '-v7.3');
rebuilt = true;
end

function [pmfPadByModel, matchPropByModel, rebuilt] = loadOrBuildAllModelPostPredSummary( ...
  summaryPath, forceRegen, ...
  storageDir, dataName, engine, ...
  cognitiveStems, contaminantStems, cognitiveJagsNames, modelLong, ...
  data, d, nP, nT, nPairs, nModels, nTpMat, nPredSamples, rngSeed, ...
  choiceThreshold)
%LOADORBUILDALLMODELPOSTPREDSUMMARY  Cache PP PMFs + forced-choice match for all models.

needBuild = forceRegen || ~isfile(summaryPath);

if ~needBuild
  S = load(summaryPath);
  needBuild = ~(isfield(S, 'pmfPadByModel') && isfield(S, 'matchPropByModel') && ...
    isfield(S, 'sourcePaths') && isfield(S, 'sourceMtimes') && ...
    isfield(S, 'nPredSamples') && isfield(S, 'rngSeed') && ...
    isfield(S, 'nP') && isfield(S, 'nModels') && isfield(S, 'choiceThreshold'));
  if ~needBuild
    needBuild = ~(S.nPredSamples == nPredSamples && S.rngSeed == rngSeed && ...
      S.nP == nP && S.nModels == nModels && S.choiceThreshold == choiceThreshold);
  end
  if ~needBuild
    curMt = nan(numel(S.sourcePaths), 1);
    for k = 1:numel(S.sourcePaths)
      curMt(k) = fileMtime(S.sourcePaths{k});
    end
    needBuild = any(~isfinite(curMt)) || ~isequal(curMt(:), S.sourceMtimes(:));
  end
  % Rebuild when a previously missing model fit has since appeared
  if ~needBuild
    for mi = 1:nModels
      if mi <= 8
        stem = cognitiveStems{mi};
      else
        stem = contaminantStems{mi - 8};
      end
      fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', stem, dataName, engine));
      if isfile(fpath) && ~any(strcmp(S.sourcePaths, fpath))
        needBuild = true;
        break;
      end
    end
  end
  if ~needBuild
    pmfPadByModel = S.pmfPadByModel;
    matchPropByModel = S.matchPropByModel;
    rebuilt = false;
    return;
  end
end

fprintf('Building all-model post-pred summary from full chains...\n');
rng(rngSeed);

pmfPadByModel = cell(nModels, 1);
matchPropByModel = nan(nModels, nP);
sourcePaths = {};
sourceMtimes = [];

for mi = 1:nModels
  if mi <= 8
    stem = cognitiveStems{mi};
  else
    stem = contaminantStems{mi - 8};
  end
  fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', stem, dataName, engine));
  pmfPadByModel{mi} = cell(1, nP);
  if ~isfile(fpath)
    warning('Missing fit for model %s (%s)', modelLong{mi}, fpath);
    continue;
  end
  fprintf('Loading %s\n', fpath);
  L = load(fpath, 'chains');
  chains = L.chains;
  sourcePaths{end + 1} = fpath; %#ok<AGROW>
  sourceMtimes(end + 1, 1) = fileMtime(fpath); %#ok<AGROW>

  for pp = 1:nP
    if mi <= 8
      Th = hierarchicalExecutionThetaDraws( ...
        chains, data, cognitiveJagsNames{mi}, pp);
      pmfPadByModel{mi}{pp} = postPredictiveCountPmfFromTheta( ...
        Th, d.pair(pp, :), nPairs, nTpMat(pp, :), nPredSamples);
    else
      pmfPadByModel{mi}{pp} = postPredictiveCountPmfFromMonitoredTheta( ...
        chains, pp, d.pair(pp, :), nT, nPairs, nTpMat(pp, :), nPredSamples);
    end

    meanTh = mapModelMeanTheta(chains, data, cognitiveJagsNames, mi, pp, nT);
    if isempty(meanTh) || numel(meanTh) ~= nT
      continue;
    end
    obs = d.LL(pp, :);
    forcedLL = nan(size(meanTh));
    forcedLL(meanTh > choiceThreshold) = 1;
    forcedLL(meanTh < choiceThreshold) = 0;
    valid = isfinite(forcedLL) & isfinite(obs);
    if any(valid)
      matchPropByModel(mi, pp) = mean(forcedLL(valid) == obs(valid));
    end
  end
  clear chains L;
end

save(summaryPath, ...
  'pmfPadByModel', 'matchPropByModel', ...
  'sourcePaths', 'sourceMtimes', ...
  'nPredSamples', 'rngSeed', 'nP', 'nModels', 'choiceThreshold', '-v7.3');
rebuilt = true;
end

function [mapModel, mapProb, matchProp, meanThetaByParticipant, rebuilt] = ...
  loadOrBuildMapForcedChoiceMatchSummary( ...
    summaryPath, forceRegen, ...
    mixPath, storageDir, dataName, engine, ...
    mixtureName, cognitiveStems, contaminantStems, cognitiveJagsNames, modelLong, ...
    data, d, nP, nT, nModels, choiceThreshold)
%LOADORBUILDMAPFORCEDCHOICEMATCHSUMMARY  Cache MAP forced-choice match rates.
%
%   MAP identity comes from mixPath. Hierarchical theta comes from
%   cognitiveStems / contaminantStems (pass MuPrecHalf/Double stems for
%   prior-robustness match rates).

needBuild = forceRegen || ~isfile(summaryPath);

if ~needBuild
  S = load(summaryPath);
  needBuild = ~(isfield(S, 'mapModel') && isfield(S, 'mapProb') && ...
    isfield(S, 'matchProp') && isfield(S, 'meanThetaByParticipant') && ...
    isfield(S, 'sourcePaths') && isfield(S, 'sourceMtimes') && ...
    isfield(S, 'choiceThreshold') && isfield(S, 'nP') && isfield(S, 'mixtureName'));
  if ~needBuild
    needBuild = ~(S.choiceThreshold == choiceThreshold && S.nP == nP && ...
      strcmp(S.mixtureName, mixtureName) && ...
      ~isempty(S.sourcePaths) && strcmp(S.sourcePaths{1}, mixPath));
  end
  if ~needBuild
    curMt = nan(numel(S.sourcePaths), 1);
    for k = 1:numel(S.sourcePaths)
      curMt(k) = fileMtime(S.sourcePaths{k});
    end
    needBuild = any(~isfinite(curMt)) || ~isequal(curMt(:), S.sourceMtimes(:));
  end
  % Rebuild when a previously missing MAP-model fit has since appeared on disk
  if ~needBuild && isfield(S, 'mapModel')
    for mi = unique(S.mapModel(:))'
      if ~(isfinite(mi) && mi >= 1 && mi <= 8)
        continue;
      end
      fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', ...
        cognitiveStems{mi}, dataName, engine));
      if isfile(fpath) && ~any(strcmp(S.sourcePaths, fpath))
        needBuild = true;
        break;
      end
    end
  end
  if ~needBuild
    mapModel = S.mapModel;
    mapProb = S.mapProb;
    matchProp = S.matchProp;
    meanThetaByParticipant = S.meanThetaByParticipant;
    rebuilt = false;
    return;
  end
end

fprintf('Building MAP forced-choice match summary from full chains...\n');

fprintf('Loading mixture %s\n', mixPath);
mixS = load(mixPath, 'chains');
postModelProb = posteriorModelProbsFromZ(mixS.chains, nP, nModels);
[mapProb, mapModel] = max(postModelProb, [], 1);
clear mixS;

uniqueMaps = unique(mapModel(~isnan(mapModel)));
chainsByModel = cell(nModels, 1);
modelPaths = {};
modelMtimes = [];
for mi = uniqueMaps(:)'
  if mi <= 8
    stem = cognitiveStems{mi};
  else
    stem = contaminantStems{mi - 8};
  end
  fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', stem, dataName, engine));
  if ~isfile(fpath)
    if mi <= 8 && (contains(stem, 'MuPrecHalf') || contains(stem, 'MuPrecDouble'))
      fprintf('Skipping MAP model %s — prior-robustness fit not ready:\n  %s\n', ...
        modelLong{mi}, fpath);
      continue;
    end
    warning('Missing fit for MAP model %s (%s)', modelLong{mi}, fpath);
    continue;
  end
  fprintf('Loading %s\n', fpath);
  L = load(fpath, 'chains');
  chainsByModel{mi} = L.chains;
  modelPaths{end + 1} = fpath; %#ok<AGROW>
  modelMtimes(end + 1, 1) = fileMtime(fpath); %#ok<AGROW>
end

matchProp = nan(1, nP);
meanThetaByParticipant = cell(1, nP);
for pp = 1:nP
  mi = mapModel(pp);
  if ~(isfinite(mi) && mi >= 1 && mi <= nModels) || isempty(chainsByModel{mi})
    continue;
  end
  % Prefer stem actually loaded (may be MuPrec*) for theta mapping
  if mi <= 8
    jagsNames = cognitiveStems;
  else
    jagsNames = cognitiveJagsNames;
  end
  meanTh = mapModelMeanTheta( ...
    chainsByModel{mi}, data, jagsNames, mi, pp, nT);
  if isempty(meanTh) || numel(meanTh) ~= nT
    continue;
  end
  meanThetaByParticipant{pp} = meanTh;
  obs = d.LL(pp, :);
  forcedLL = nan(size(meanTh));
  forcedLL(meanTh > choiceThreshold) = 1;
  forcedLL(meanTh < choiceThreshold) = 0;
  valid = isfinite(forcedLL) & isfinite(obs);
  if any(valid)
    matchProp(pp) = mean(forcedLL(valid) == obs(valid));
  end
end
clear chainsByModel;

sourcePaths = [{mixPath}; modelPaths(:)];
sourceMtimes = [fileMtime(mixPath); modelMtimes(:)];
save(summaryPath, ...
  'mapModel', 'mapProb', 'matchProp', 'meanThetaByParticipant', ...
  'sourcePaths', 'sourceMtimes', ...
  'choiceThreshold', 'nP', 'mixtureName', '-v7.3');
rebuilt = true;
end

function meanTh = mapModelMeanTheta(chains, data, cognitiveJagsNames, mi, pp, nT)
%MAPMODELMEANTHETA  Posterior-mean P(LL) per trial under MAP model mi.
if mi <= 8
  Th = hierarchicalExecutionThetaDraws(chains, data, cognitiveJagsNames{mi}, pp);
  if isempty(Th)
    meanTh = [];
    return;
  end
  meanTh = mean(Th, 1);
  return;
end

meanTh = nan(1, nT);
for j = 1:nT
  fn = sprintf('theta_%d_%d', pp, j);
  if ~isfield(chains, fn)
    meanTh = [];
    return;
  end
  meanTh(j) = mean(chains.(fn)(:));
end
end

function mt = fileMtime(path)
d = dir(path);
if isempty(d)
  mt = NaN;
else
  mt = d.datenum;
end
end

function P = posteriorModelProbsFromZ(chains, nParticipants, nModels)
P = nan(nModels, nParticipants);
zMats = codaIndexedMatrices(chains, 'z');
for pp = 1:min(nParticipants, numel(zMats))
  z = round(zMats{pp}(:));
  z = z(isfinite(z) & z >= 1 & z <= nModels);
  if isempty(z)
    continue;
  end
  cnt = histcounts(z, 0.5:(nModels + 0.5));
  P(:, pp) = cnt(:) / sum(cnt);
end
end

function fam = mapModelToUseFamily(mapModel)
%MAPMODELTOUSEFAMILY  Collapse tr/ut (indices 6/7) to one model-use family.
fam = mapModel;
fam(mapModel == 7) = 6; % unified tradeoff -> tradeoff family
end

function Pfam = combinedFamilyPosteriors(P)
%COMBINEDFAMILYPOSTERIORS  Merge Tr+UT (indices 6+7) into one family mass.
%
%   P is nModels-by-nParticipants. Returns 10-by-nP:
%   Ex, Hc, Hd, PD, DD, Tr+UT, IT, Gu, LL, SS.
Pfam = [P(1:5, :); P(6, :) + P(7, :); P(8:11, :)];
end

function tf = mixtureHasIndexedParam(chains, paramName)
tf = false;
if isempty(chains) || ~isstruct(chains)
  return;
end
try
  mats = codaIndexedMatrices(chains, paramName);
  tf = ~isempty(mats);
catch
  tf = false;
end
end

function [mu, qLo, qHi] = indexedParamSummary(chains, paramName, pp, iqrProbs)
mu = NaN;
qLo = NaN;
qHi = NaN;
if isempty(chains) || ~isstruct(chains)
  return;
end
try
  mats = codaIndexedMatrices(chains, paramName);
catch
  return;
end
if pp < 1 || pp > numel(mats) || isempty(mats{pp})
  return;
end
s = mats{pp}(:);
s = s(isfinite(s));
if isempty(s)
  return;
end
mu = mean(s);
qs = quantile(s, iqrProbs);
qLo = qs(1);
qHi = qs(2);
end

function mu = indexedParamMean(chains, paramName, pp)
[mu, ~, ~] = indexedParamSummary(chains, paramName, pp, [0.25 0.75]);
end

function [mu, qLo, qHi] = paramSummaryFromMapModel( ...
  indivChains, mi, indivField, pp, modelParamSpecs, iqrProbs)
mu = NaN;
qLo = NaN;
qHi = NaN;
if ~(isfinite(mi) && mi >= 1 && mi <= numel(modelParamSpecs))
  return;
end
if mi > numel(indivChains) || isempty(indivChains{mi})
  return;
end
specs = modelParamSpecs{mi};
hasField = false;
for s = 1:numel(specs)
  if strcmp(specs{s}{2}, indivField)
    hasField = true;
    break;
  end
end
if ~hasField
  return;
end
[mu, qLo, qHi] = indexedParamSummary(indivChains{mi}, indivField, pp, iqrProbs);
end

function mu = paramMeanFromMapModel(indivChains, mi, indivField, pp, modelParamSpecs)
[mu, ~, ~] = paramSummaryFromMapModel( ...
  indivChains, mi, indivField, pp, modelParamSpecs, [0.25 0.75]);
end

function pmfPad = postPredictiveCountPmfFromMonitoredTheta(chains, pp, pairRow, nT, nPairs, nTpRow, nPredSamples)
% Contaminant models (guess / LL / SS) store monitored theta_pp_j.
nkMax = max(nTpRow);
if nkMax <= 0
  pmfPad = [];
  return;
end

pmfPad = nan(nPairs, nkMax + 1);
for j = 1:nT
  fn = sprintf('theta_%d_%d', pp, j);
  if ~isfield(chains, fn)
    pmfPad = [];
    return;
  end
end

fn1 = sprintf('theta_%d_%d', pp, 1);
nS = numel(chains.(fn1)(:));
Th = nan(nS, nT);
for j = 1:nT
  fn = sprintf('theta_%d_%d', pp, j);
  Th(:, j) = chains.(fn)(:);
end

if nS > nPredSamples && nPredSamples > 0
  idx = randperm(nS, nPredSamples);
  Th = Th(idx, :);
end

Th = min(max(Th, 0), 1);
useS = size(Th, 1);
cntSim = nan(useS, nPairs);
for s = 1:useS
  ySim = rand(1, nT) < Th(s, :);
  for p = 1:nPairs
    m = pairRow == p;
    if any(m)
      cntSim(s, p) = sum(ySim(m));
    end
  end
end

for p = 1:nPairs
  nk = nTpRow(p);
  if nk <= 0
    continue;
  end
  col = cntSim(:, p);
  col = col(isfinite(col));
  if isempty(col)
    continue;
  end
  for k = 0:nk
    pmfPad(p, k + 1) = mean(col == k);
  end
end
end

function stem = cognitiveMuPrecStem(baseStem, precLabel)
% Insert MuPrecHalf / MuPrecDouble before _entrop.
stem = regexprep(baseStem, '_entrop$', ['MuPrec' precLabel '_entrop']);
if strcmp(stem, baseStem)
  error('cognitiveMuPrecStem:badStem', ...
    'Expected stem ending in _entrop, got: %s', baseStem);
end
end

function tf = hierarchicalStemAvailable(storageDir, stem, dataName, engine)
tf = isfile(fullfile(storageDir, sprintf('%s_%s_%s.mat', stem, dataName, engine)));
end

function chains = loadHierarchicalChains(storageDir, stem, dataName, engine)
fpath = fullfile(storageDir, sprintf('%s_%s_%s.mat', stem, dataName, engine));
if ~isfile(fpath)
  error(['Missing hierarchical fit %s\n' ...
    'Run models/runHierarchicalExecutionPriorRobustness.m for half/double fits.'], ...
    fpath);
end
fprintf('Loading %s\n', fpath);
L = load(fpath, 'chains');
chains = L.chains;
end

function writeModelParticipantTable(outTex, values, mapModel, modelShort, caption, label, asPercent)
%WRITEMODELPARTICIPANTTABLE  APA sidewaystable: participants x models.
nModels = numel(modelShort);
nP = size(values, 2);
fid = fopen(outTex, 'w');
if fid < 0
  error('Could not write %s', outTex);
end
cleanup = onCleanup(@() fclose(fid));

fprintf(fid, '%% Auto-generated by drawFiguresEntrop (modelAgreementTables) — do not edit by hand.\n');
fprintf(fid, '\\begin{sidewaystable}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{%s}\n', caption);
fprintf(fid, '\\label{%s}\n', label);
fprintf(fid, '\\begin{adjustbox}{max width=\\textheight}\n');
fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('r', 1, nModels));
fprintf(fid, '\\toprule\n');
fprintf(fid, 'Participant');
for mi = 1:nModels
  fprintf(fid, ' & %s', modelShort{mi});
end
fprintf(fid, ' \\\\\n\\midrule\n');

for pp = 1:nP
  fprintf(fid, '%s', char('A' + pp - 1));
  miMap = mapModel(pp);
  for mi = 1:nModels
    v = values(mi, pp);
    if ~isfinite(v)
      cellStr = '---';
    elseif asPercent
      cellStr = sprintf('%.0f', round(100 * v));
    else
      cellStr = sprintf('%.2f', v);
    end
    if isfinite(miMap) && mi == miMap
      fprintf(fid, ' & \\textbf{%s}', cellStr);
    else
      fprintf(fid, ' & %s', cellStr);
    end
  end
  fprintf(fid, ' \\\\\n');
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{adjustbox}\n');
fprintf(fid, '\\end{sidewaystable}\n');
fprintf('Wrote %s\n', outTex);
end

function writeModelParticipantCsv(csvPath, values, modelShort)
%WRITEMODELPARTICIPANTCSV  Participant x model values as CSV.
nModels = numel(modelShort);
nP = size(values, 2);
hdr = [{'Participant'}, modelShort];
rows = cell(nP, nModels + 1);
for pp = 1:nP
  rows{pp, 1} = char('A' + pp - 1);
  for mi = 1:nModels
    rows{pp, mi + 1} = values(mi, pp);
  end
end
writecell([hdr; rows], csvPath);
fprintf('Wrote %s\n', csvPath);
end
