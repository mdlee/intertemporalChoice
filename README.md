# Intertemporal choice: evaluating integration rules

Code, data, and figures for an individual-level Bayesian comparison of eight intertemporal-choice integration rules (exponential, hyperbolic, hyperboloid, proportional differences, direct differences, tradeoff, unified tradeoff, and ITCH), plus three contaminant models, using Grünwald entropification as the decision rule.

The paper figures live under `results/`. MCMC chain files are not in the repository (they are large); re-running the MATLAB scripts below regenerates them under `models/storage/`.

## Layout

| Path | What it is |
|------|------------|
| `data/` | Trial-level CSVs (`sujeto_*_tiempo.csv`) and `intertemporalChoice.mat` |
| `models/jags/` | JAGS implementations used in the paper |
| `models/` | MATLAB drivers for hierarchical fits, the latent mixture, and prior robustness |
| `general/` | Small Trinity/JAGS and figure helpers (R-hat, coda, layout) |
| `results/` | Publication figures and the extra panels linked from the paper |
| `results/drawFiguresEntrop.m` | Builds every figure in `results/` from stored fits |

## Paper links in this repository

- [JAGS models](https://github.com/mdlee/intertemporalChoice/tree/main/models/jags)
- [Latent-mixture posteriors, including half/double prior versions](https://github.com/mdlee/intertemporalChoice/tree/main/results/latentMixture)
- [Descriptive adequacy and posterior-predictive robustness](https://github.com/mdlee/intertemporalChoice/tree/main/results/descriptiveAdequacy)
- [Parameter inferences and prior robustness](https://github.com/mdlee/intertemporalChoice/tree/main/results/parameterInferences)

## Requirements

- MATLAB
- [Trinity](https://github.com/jooh/trinity) on the MATLAB path (`callbayes` / JAGS)
- JAGS

From a MATLAB session whose current folder is this repository:

```matlab
addpath(fullfile(pwd, 'general'));
addpath(fullfile(pwd, 'models'));
```

`data/intertemporalChoice.mat` is already prepared. To rebuild it from the CSVs, run `data/parse_1.m`.

## Reproducing the analyses

Fits write `models/storage/` (gitignored). Half/double prior JAGS files are generated into `models/jagsRobust/` (also gitignored) by `generateHierarchicalMuPrecJags.m`.

1. **Hierarchical cognitive models** (and half/double prior-robustness variants):

   ```matlab
   cd models
   runHierarchicalExecutionSequential
   ```

   Or fit one model at a time with `runHierarchicalExecutionModel`.

2. **Latent mixture** (11 components: eight cognitive models + guess / larger-later / smaller-sooner):

   ```matlab
   runLatentMixtureSequential
   ```

3. **Mixture prior robustness** (half and double group-mean prior precision):

   ```matlab
   runLatentMixtureRobustnessSequential
   ```

4. **Figures** (after the `.mat` fits exist):

   ```matlab
   cd ../results
   drawFiguresEntrop
   ```

   Toggle which panels to rebuild in `analysisList` near the top of `drawFiguresEntrop.m`.

`runHierarchicalExecutionPriorRobustness` / `runHierarchicalExecutionPriorRobustnessMAP` refit the separate hierarchical models under half/double priors if you are not using the sequential driver.

## Results folders

`results/basic/` has the experimental design and choice counts. `results/latentMixture/` has posterior model probabilities with descriptive-adequacy markers (original, half, and double priors). `results/descriptiveAdequacy/` has MAP posterior predictives, the misfit gallery, and forced-choice robustness. `results/parameterInferences/` has selected-model parameters, all-participant parameters under each model, and prior-robustness scatters.

## License / citation

Please cite the paper if you use the data or code. See `general/trinity_license.txt` for Trinity components vendored in `general/`.
