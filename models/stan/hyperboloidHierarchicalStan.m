%% hyperboloidHierarchicalStan — hierarchical model via Stan (CmdStan / Trinity)
% Same specification as ../hyperboloidHierarchical.m; chains under stan/storage/.
cd(fileparts(mfilename('fullpath')));
runHierarchicalModelStan('hyperboloidHierarchical', {'kappa', 'tau', 'sigma'});
