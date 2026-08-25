%% hyperbolicHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../hyperbolicHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('hyperbolicHierarchical', {'kappa', 'sigma'});
