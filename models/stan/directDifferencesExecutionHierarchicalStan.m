%% directDifferencesExecutionHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../directDifferencesExecutionHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('directDifferencesExecutionHierarchical', {'delta', 'omega', 'epsilon'});
