%% hyperbolicExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../hyperbolicExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('hyperbolicExecutionHierarchical', {'kappa', 'epsilon'});
