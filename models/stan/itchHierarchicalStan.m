%% itchHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../itchHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('itchHierarchical', {'beta0', 'betaRA', 'betaRR', 'betaTA', 'betaTR', 'sigma'});
