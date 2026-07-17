# Local general utilities

Local copies of shared MATLAB helpers used by the intertemporal-choice analyses,
so the project does not depend on a machine-wide `MatlabAddOns` path for these.

## Contents

| File / folder | Role |
|---|---|
| `findKeepChains.m` | Subset MCMC chains to meet an R-hat target |
| `gelmanrubin.m` | Gelman–Rubin R-hat (Trinity) |
| `get_matrix_from_coda.m`, `grtable.m`, `codatable.m` | Trinity coda summaries |
| `+trinity/` | Minimal Trinity package pieces (`select_fields`, `error_tag`) |
| `pantoneColors.mat` | Color palette used by figure scripts |
| `setFigure.m`, `subplotArrange.m`, `moveAxis.m`, `suplabel.m` | Figure layout helpers |

`callbayes` / full Trinity JAGS engine is **not** vendored here — keep Trinity on the MATLAB path for sampling.

Sources: `~/Dropbox/MatlabAddOns` and `trinity-master` (see `trinity_license.txt`).
