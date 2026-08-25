%% latentMixtureHierarchicalPrecision_1plusStan — latent mixture via Stan
% probit (1plus). Discrete z is marginalized; see latentMixtureHierarchicalPrecision_1plus.stan generated quantities.
clear; close all;
cd(fileparts(mfilename('fullpath')));
addpath(fileparts(fileparts(mfilename('fullpath'))));
modelName = 'latentMixtureHierarchicalPrecision_1plus';
[~, d] = prepareIntertemporalChoiceData();
initGenerator = @() struct('sigmaHY', max(randn(d.nParticipants, 1) + 1, 0.5));
runHierarchicalLatentMixtureModelStan(modelName, initGenerator);
