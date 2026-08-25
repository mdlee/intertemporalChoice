%% exponentialExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../exponentialExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('exponentialExecutionHierarchical', {'kappa', 'epsilon'});
