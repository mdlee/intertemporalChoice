# Stan hierarchical models (MATLAB drivers)

Stan reimplementation of the hierarchical execution, probit, and latent-mixture models in the parent `models/` folder. **MATLAB still orchestrates** data prep, MCMC loops, diagnostics, and figures via Trinity `callbayes('stan', ...)`.

## Layout

| Path | Role |
|------|------|
| `*.stan` | Source models (`#include "common_functions.stan"` for editing only) |
| `common_functions.stan` | Shared decision rules; merged into staged copies under `tmp/` before compile |
| `runHierarchicalModelStan.m` | Single-model runner (mirrors `runHierarchicalExecutionModel.m`) |
| `runHierarchicalLatentMixtureModelStan.m` | Latent mixture runner (mirrors `runHierarchicalLatentMixtureModel.m`) |
| `*Stan.m` | Thin entry scripts (one per model) |
| `runAllHierarchicalExecutionModelsStan.m` | Batch fit all 16 single hierarchical models |
| `storage/` | `{modelName}_intertemporalChoice_stan.mat` |
| `figures/` | `{modelName}_thin{N}.png`, `_final.png`, etc. |
| `tmp/{modelName}/` | CmdStan working directory per fit |

Parent-folder utilities (`prepareIntertemporalChoiceData`, `hierarchicalFigurePath`, `plotParticipantParameterCIs`, `hierarchicalExecutionThetaDraws`, MCMC state helpers) are on the path via `addpath` in the runners.

## Requirements

- MATLAB with Trinity (`callbayes`) and CmdStan (e.g. `stan_main_dir` in Trinity preferences).
- Intertemporal-choice data as for the JAGS pipeline.

## Usage

```matlab
cd('/path/to/models/stan');
run('exponentialExecutionHierarchicalStan');   % one model
run('runAllHierarchicalExecutionModelsStan'); % all 16 singles
run('latentMixtureHierarchicalPrecision_2plusStan');
```

Options match the JAGS runners (`preLoad`, `resetThin`, `nChains`, `rhatCritical`, …). Trinity is called with `allowunderscores` because Stan hyperparameters use names like `mu_kappa` (Trinity otherwise reserves `_` for indexed CODA names).

## Latent mixture (manuscript 11-component)

The manuscript analysis is the 11-component model in
`../jags/latentMixtureHierarchicalPrecision_entrop_jags.txt`:

- one discrete `z[i]` per participant (not per trial)
- 8 hierarchical cognitive models plus Guess / LL / SS contaminants
- Grünwald entropification: `P(LL) = inv_logit(±w)` after a hard `step()` comparison

Stan counterpart: `latentMixtureHierarchicalPrecision_entrop.stan`. Discrete `z` is
**marginalized** with a participant-level `log_sum_exp` over the 11 component
likelihoods. Generated quantities draw `z[i] ~ categorical_logit(lp)` (a posterior
assignment draw, not an argmax).

μ-precision robustness uses the same `.stan` file with `data.muPrecScale` = 0.5 or 2.
Sequential driver (from `models/`):

```matlab
cd('/path/to/models');
run('runLatentMixtureSequentialStan.m');  % base, half, double
```

Start that **after** `runHierarchicalExecutionSequentialStan` finishes (warm-start
from the 8 cognitive Stan fits). A watcher does this:

```bash
models/stan/logs/queueLatentMixtureAfterHierStan.sh
```

Gate: best subset of at least 4 of 6 chains, max participant z R-hat ≤ 1.10, n ≥ 10,000.

Legacy `latentMixtureHierarchicalPrecision_1plus.stan` / `_2plus.stan` are earlier
ports (probit/ε noise, per-trial mixture) and are **not** the manuscript model.

## Priors

Single hierarchical probit `.stan` files use the same structural hyperpriors as execution (e.g. `mu_sigma ~ normal(0.05, 0.1)`). Latent-mixture Stan files follow the 2plus/1plus JAGS structure with bounded parameters; noise hyperpriors use the 0.05-scale pattern used in the 2plus Stan port.
