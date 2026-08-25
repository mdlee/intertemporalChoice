%% proportionalDifferencesHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../proportionalDifferencesHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('proportionalDifferencesHierarchical', {'delta', 'sigma'});
