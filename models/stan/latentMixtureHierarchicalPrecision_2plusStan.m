%% latentMixtureHierarchicalPrecision_2plusStan — latent mixture via Stan
% execution (2plus). Discrete z is marginalized; see latentMixtureHierarchicalPrecision_2plus.stan generated quantities.
clear; close all;
cd(fileparts(mfilename('fullpath')));
addpath(fileparts(fileparts(mfilename('fullpath'))));
modelName = 'latentMixtureHierarchicalPrecision_2plus';
[~, d] = prepareIntertemporalChoiceData();
initGenerator = @() struct('epsilonHY', rand(d.nParticipants, 1) * 0.4 + 0.05);
runHierarchicalLatentMixtureModelStan(modelName, initGenerator);
