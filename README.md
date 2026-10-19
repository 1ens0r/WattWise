# WattWise

Household energy prediction with linear models, gradient descent and statistical
inference, applied to the UCI **Appliances Energy Prediction** dataset in R.

The project predicts appliance energy consumption (watt-hours) 10 minutes ahead
from recent consumption, indoor climate, outdoor weather and calendar time, and
scores every model on a **chronological holdout** period using the same protocol.
All models are linear (least squares): penalised, expanded, lag-structured and
partitioned variants of the same underlying estimator.

Every number in `results/`, in `results.md` and in `report.md` is produced by the
scripts in this repository.

---

## Project layout

```text
WattWise_Phase1/
├── run_all.R                       # run the whole analysis end to end
├── results.md                      # generated head-to-head benchmark write-up
├── report.md                       # generated project report (the full narrative)
├── LICENSE                         # MIT
├── eda/                            # data preparation and exploratory analysis
│   └── 00_prepare_data.R
├── benchmarks/                     # baselines, model harness, solver checks
│   ├── 01_naive_baselines.R
│   ├── 02_model_benchmark.R
│   └── 03_solver_checks.R
├── models/                         # the six benchmark models
│   ├── 00_shared.R                 # shared features, blocked CV, solvers
│   ├── 01_ridge.R
│   ├── 02_elastic_net.R
│   ├── 03_lasso.R
│   ├── 04_polynomial_ridge.R
│   ├── 05_distributed_lag.R
│   └── 06_phmr.R
├── archive/models/                 # superseded phase-1 models (still runnable)
├── results/
│   ├── md/                         # written by the scripts; one summary per study
│   ├── csv/                        # machine-readable tables written by the scripts
│   └── plots/                      # figures (populated by the EDA work)
├── data/
│   ├── raw/                        # downloaded source CSV (read-only input)
│   └── processed/                  # cleaned data and the chronological split
└── WattWise_Phase1.Rproj
```

### Conventions

- Folder = stage: `eda/` prepares and describes data, `benchmarks/` holds
  baselines, the benchmark harness and verification, `models/` holds one file per
  model.
- Script names are `NN_snake_case.R`; the numeric prefix is the run order.
- Scripts write their own summaries to `results/md/` and their tables to
  `results/csv/`. Documentation never restates numbers by hand, so there is
  nothing to keep in sync.
- Every script resolves the project root itself and runs from any working
  directory.
- The models are implemented in **base R only**: ridge, lasso and elastic net are
  solved from scratch, and `benchmarks/03_solver_checks.R` verifies them against
  independent optimality conditions.
- Tuning is the parallel part of the pipeline. Each model is tuned by
  cross-validation, and those folds are independent, so the harness spreads them
  over the machine's cores (up to 8) with a PSOCK cluster, and falls back to a
  single process automatically if a cluster cannot be created. The score for a
  model never depends on the number of cores: the folds are combined by summing
  in a fixed order. Set `options(wattwise.cores = 1L)` to force serial execution.

---

## The six models

| # | Model | Why it suits this dataset | Reference |
| - | --- | --- | --- |
| 1 | Ridge | Indoor and outdoor sensors are strongly correlated, which makes plain OLS coefficients unstable; L2 shrinkage trades bias for variance. | Hoerl & Kennard (1970), *Technometrics* 12(1), 55–67 |
| 2 | Elastic net | The nine indoor temperature and nine humidity channels are collinear: the L1+L2 mix keeps correlated groups together instead of picking one arbitrarily. | Zou & Hastie (2005), *JRSS-B* 67(2), 301–320 |
| 3 | Lasso | Acts as an embedded feature selector, giving an explicit, defensible shortlist of the 28 sensor and calendar channels. | Tibshirani (1996), *JRSS-B* 58(1), 267–288 |
| 4 | Polynomial ridge | Appliance use responds to heating/cooling thresholds and humidity–temperature interactions; a degree-2 expansion adds curvature while staying linear in its coefficients. | Hastie, Tibshirani & Friedman (2009), *ESL* 2nd ed., Ch. 3.1/3.4 |
| 5 | Distributed lag (OLS) | Thermal inertia means weather effects are not instantaneous; lagged outdoor channels estimate the whole lag profile and test whether history adds information. | Almon (1965), *Econometrica* 33(1), 178–196 |
| 6 | PHMR | One global linear relationship is a strong assumption for a house that behaves differently through the day; PHMR fits related per-regime regressions with partial pooling. | Ning & Yoon (2024), *Indoor Air* 2024(1), 5595459 |

Full citations, including the solver reference (Friedman, Hastie & Tibshirani
2010) and the dataset-specific published benchmark (Lu, Lu & Kassymova 2026), are
in the references section of `results.md` and at the top of each model file.

## Benchmark protocol

Identical for every model, and asserted in code:

1. **Target.** `Appliances` at *t+1*, one 10-minute step ahead.
2. **Information set.** Features come from *t* or earlier only. `lights` is
   excluded from the models, and the dataset's random columns `rv1`/`rv2` are
   never used.
3. **Split.** The chronological 80/20 holdout from `eda/00_prepare_data.R`. Every
   training timestamp precedes every holdout timestamp.
4. **No holdout leakage.** Standardisation statistics and every hyperparameter
   come from training rows only.
5. **Time-aware tuning.** Cross-validation uses contiguous forward-chaining
   (blocked) folds, never random splits, because random splits leak across time
   in an autocorrelated series.
6. **Determinism.** Closed-form solves and zero-initialised coordinate descent, so
   results reproduce exactly with no random seed.

## Requirements

- R 4.1 or newer. Verified with R 4.6.1 on Windows.
- **Base R only** — no packages to install.
- Internet access only if `data/raw/energydata_complete.csv` is absent (the prep
  script downloads it from UCI automatically).
- The pipeline uses the available cores for tuning (see Conventions). On a
  single-core machine it still runs, just more slowly.

## How to run

Open `WattWise_Phase1.Rproj` in RStudio (which sets the working directory to the
project root), then either run everything:

```r
source("run_all.R")
```

or a single stage:

```r
source("eda/00_prepare_data.R")            # data -> data/processed + quality report
source("benchmarks/01_naive_baselines.R")  # naive reference scores
source("benchmarks/03_solver_checks.R")    # verify the from-scratch solvers
source("benchmarks/02_model_benchmark.R")  # six models head to head -> results.md
```

From a terminal, from the project root:

```sh
Rscript run_all.R
```

`run_all.R` also rebuilds `results/md/key_results.md` and `report.md` by reading
the generated tables back from `results/csv/`.

Tuning is CPU-bound and parallel across folds. A full end-to-end run takes about
two minutes on eight logical cores (most of it the elastic-net path).

## Scripts at a glance

| Script | Purpose | Writes |
| --- | --- | --- |
| `run_all.R` | Runs all stages in order and rebuilds the headline index and the report. | `results/md/key_results.md`, `report.md` |
| `eda/00_prepare_data.R` | Downloads, validates, cleans and chronologically splits the data; records the data-quality report. | `data/processed/*.csv`, `results/md/data_quality.md` |
| `benchmarks/01_naive_baselines.R` | Persistence (`t-1`) and seasonal naive (previous day) reference forecasts. | `results/csv/naive_baselines.csv`, `results/md/naive_baselines.md` |
| `benchmarks/03_solver_checks.R` | Verifies the from-scratch solvers: OLS normal equations, ridge against an independent QR solve, lasso/elastic-net KKT conditions, PHMR block optimality, determinism. | `results/md/solver_checks.md`, `results/csv/solver_checks.csv` |
| `benchmarks/02_model_benchmark.R` | Trains all six models plus three references on the identical holdout, tunes on training data only, and writes up the comparison. | `results.md`, `results/csv/model_*.csv` |
| `models/00_shared.R` | Shared library: leakage-free feature engineering, blocked CV, solvers, metrics. Sourced by the harness. | – |
| `models/01_ridge.R` … `models/06_phmr.R` | One `wattwise_model` definition per file (tune + fit), consumed by the harness. | – |

## Results

Start with **`report.md`** — the written project report: problem statement, data,
methods with citations, the results table, discussion and limitations. The
mechanical detail lives in **`results.md`** — the generated head-to-head benchmark
(protocol, per-model rationale, paired comparisons, stability, limitations).

| File | Contents |
| --- | --- |
| `report.md` | The project report: narrative, methods, results table, discussion, references. |
| `results.md` | The six-model benchmark write-up, ranked, with the winner named. |
| `results/md/key_results.md` | Index: split summary, benchmark table, naive baselines. |
| `results/md/data_quality.md` | Dimensions, schema check, cleaning counts, split boundary. |
| `results/md/naive_baselines.md` | Persistence vs seasonal naive, with the benchmark to beat. |
| `results/md/solver_checks.md` | Independent verification of the from-scratch solvers. |
| `results/csv/solver_checks.csv` | The same verification as a table. |
| `results/csv/model_benchmark.csv` | One row per model: RMSE, MAE, R², tuning, runtime. |
| `results/csv/model_predictions.csv` | Per-timestamp predictions, every model. |
| `results/csv/model_coefficients.csv` | Standardised coefficients (top 60 for the 464-term model). |
| `results/csv/model_paired_tests.csv` | Paired Wilcoxon comparisons against the best model. |
| `results/csv/model_stability.csv` | RMSE by calendar month for the leading models. |

## Research references

The six methods were chosen for this dataset rather than for their citation
count. Each is implemented from the source below, and the reason it suits this
data is in the model file header and in `results.md`.

| # | Method | Reference and link |
| - | --- | --- |
| 1 | Ridge regression | Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. *Technometrics*, 12(1), 55–67. <https://doi.org/10.1080/00401706.1970.10488634> |
| 2 | Elastic net | Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. *JRSS-B*, 67(2), 301–320. <https://doi.org/10.1111/j.1467-9868.2005.00503.x> |
| 3 | Lasso | Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. *JRSS-B*, 58(1), 267–288. <https://doi.org/10.1111/j.2517-6161.1996.tb02080.x> |
| 4 | Polynomial regression with regularisation | Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The Elements of Statistical Learning* (2nd ed.), Ch. 3.1 and 3.4. Springer. <https://hastie.su.domains/ElemStatLearn/> · Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*, Sec. 3.1. Springer. <https://doi.org/10.1007/978-0-387-45528-0> |
| 5 | Distributed-lag regression | Almon, S. (1965). The distributed lag between capital appropriations and expenditures. *Econometrica*, 33(1), 178–196. <https://doi.org/10.2307/1911894> |
| 6 | PHMR | Ning, S., & Yoon, H. (2024). A New Model for Building Energy Modeling and Management Using Predictive Analytics: Partitioned Hierarchical Multitask Regression (PHMR). *Indoor Air*, 2024(1), 5595459. <https://doi.org/10.1155/2024/5595459> |

Supporting references:

- **Coordinate descent solver** for lasso and elastic net: Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software*, 33(1), 1–22. <https://doi.org/10.18637/jss.v033.i01>
- **Dataset-specific published benchmark** (used for comparison, not for tuning): Lu, H., Lu, H., & Kassymova, G. K. (2026). Leakage-Aware Benchmarking of Machine-Learning Models for Residential Appliance Energy Prediction Using a Public Smart-Home Dataset. *Progress in Energy and Environment*, 32(2), 36–51. <https://doi.org/10.37934/progee.32.2.3651>
- **Related energy time-series work** (hierarchical forecast structure with a sparse penalty, the motivation behind the partitioned/hierarchical family): Pang, Y., Zhou, X., Zhang, J., Sun, Q., & Zheng, J. (2022). Hierarchical electricity time series prediction with cluster analysis and sparse penalty. *Pattern Recognition*, 126, 108555. <https://doi.org/10.1016/j.patcog.2022.108555>
- **Dataset paper**: Candanedo, L. M., Feldheim, V., & Deramaix, D. (2017). Data driven prediction models of energy use of appliances in a low-energy house. *Energy and Buildings*, 140, 81–97. <https://doi.org/10.1016/j.enbuild.2017.01.083>

## Data

- Source: UCI Machine Learning Repository, dataset 374, `energydata_complete.csv`
  (19,735 observations at 10-minute intervals, one household).
- `data/raw/` holds the downloaded file; `data/processed/` holds the cleaned
  series and the chronological 80/20 train/test split.
- The split is by time, never shuffled: every training timestamp precedes every
  test timestamp, and the scripts assert this.
- `rv1` and `rv2` are random columns shipped with the dataset and are never used
  as predictors. `lights` is excluded from the models so that predictions
  describe appliance consumption rather than total household consumption.
- Cleaning is conservative: invalid or duplicate timestamps fail validation;
  rows with missing or non-finite values are excluded; nothing is imputed or
  winsorised. The current data needs no exclusions (see the data-quality report).

## Method notes and limitations

- The evaluation is a single chronological holdout from one household over about
  five months, so scores are not estimates for other homes or seasons.
- Energy use is spiky and heavily skewed. Spanning 10 minutes, the persistence
  benchmark is strong, and a small number of switching events dominates squared
  error.
- Because a lagged consumption term is legitimately available 10 minutes ahead,
  these scores are **not** comparable with published benchmarks that omit lagged
  targets. `results.md` makes that comparison explicitly rather than implying our
  models beat the literature.
- Errors are not independent across neighbouring timestamps, so classical i.i.d.
  inference (t-tests, F-tests, ANOVA) needs care; the benchmark's paired tests are
  labelled indicative for exactly this reason.

## Status and next steps

Done: data preparation and split, naive benchmarks, from-scratch solver
verification, six linear models implemented and benchmarked head to head, and a
generated `results.md` plus the project `report.md`.

Next, all still linear-model work:

1. `eda/` descriptive statistics and distribution checks (target and residuals),
   writing plots to `results/plots/`.
2. Coefficient tables for the winning model: estimates, standard errors,
   confidence intervals, the overall F-test, R² and residual diagnostics.
3. `models/` gradient descent: batch and mini-batch descent on the same
   least-squares objective from a cold start, with loss curves compared against
   the closed-form solution.
4. Hypothesis tests with stated null and alternative hypotheses at a fixed
   significance level, including assumption checks.
5. A slide deck assembled from `report.md`, `results.md` and `results/md/`.

## Licence

This project is released under the [MIT Licence](LICENSE). Copyright (c) 2026
1ens0r.

The dataset is **not** covered by the MIT licence: it is distributed by UCI under
CC BY 4.0, so the attribution below must be kept.

## Attribution

Candanedo, L. (2017). *Appliances Energy Prediction* [Dataset]. UCI Machine
Learning Repository. DOI: <https://doi.org/10.24432/C5VC8G>. Dataset licence:
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Retain attribution
when redistributing the data or derived tables.

Related paper: Candanedo, L. M., Feldheim, V., & Deramaix, D. (2017). Data driven
prediction models of energy use of appliances in a low-energy house. *Energy and
Buildings*, 140, 81–97. <https://doi.org/10.1016/j.enbuild.2017.01.083>
