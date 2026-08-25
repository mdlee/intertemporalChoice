%% exponentialHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../exponentialHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('exponentialHierarchical', {'kappa', 'sigma'});
