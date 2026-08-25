%% proportionalDifferencesExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../proportionalDifferencesExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('proportionalDifferencesExecutionHierarchical', {'delta', 'epsilon'});
