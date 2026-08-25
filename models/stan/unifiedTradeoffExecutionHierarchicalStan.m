%% unifiedTradeoffExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../unifiedTradeoffExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('unifiedTradeoffExecutionHierarchical', {'gamma', 'tau', 'kappa', 'vartheta', 'eta', 'epsilon'});
