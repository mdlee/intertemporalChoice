%% directDifferencesHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../directDifferencesHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('directDifferencesHierarchical', {'delta', 'omega', 'sigma'});
