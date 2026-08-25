%% hyperboloidExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../hyperboloidExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('hyperboloidExecutionHierarchical', {'kappa', 'tau', 'epsilon'});
