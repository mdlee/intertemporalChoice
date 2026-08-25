%% unifiedTradeoffHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../unifiedTradeoffHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('unifiedTradeoffHierarchical', {'gamma', 'tau', 'kappa', 'vartheta', 'eta', 'sigma'});
