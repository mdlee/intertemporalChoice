%% itchExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../itchExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('itchExecutionHierarchical', {'beta0', 'betaRA', 'betaRR', 'betaTA', 'betaTR', 'epsilon'});
