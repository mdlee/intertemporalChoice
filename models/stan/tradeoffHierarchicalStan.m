%% tradeoffHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../tradeoffHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('tradeoffHierarchical', {'gamma', 'tau', 'kappa', 'vartheta', 'sigma'});
