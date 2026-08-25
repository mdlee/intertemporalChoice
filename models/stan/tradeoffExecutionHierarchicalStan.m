%% tradeoffExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../tradeoffExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('tradeoffExecutionHierarchical', {'gamma', 'tau', 'kappa', 'vartheta', 'epsilon'});
