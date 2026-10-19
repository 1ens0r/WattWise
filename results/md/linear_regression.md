# Linear regression (weather + calendar time + lag-1) -- WattWise

Generated: 2026-10-19 08:10:29 UTC
Source script: `models/01_linear_regression.R`

## Setup
- Target: `Appliances` at time *t+1* (watt-hours).
- Predictors: measurements at time *t*: `previous_appliances`, `T_out`, `RH_out`, `Windspeed`, `Visibility`, `Press_mm_hg`, `hour_sin`, `hour_cos`, `weekday_sin`, `weekday_cos`.
- Estimator: ordinary least squares via `lm()` on the training partition only; the test partition is never used for fitting.
- Evaluation period: 2016-04-30 08:30:00 to 2016-05-27 18:00:00
- Evaluated rows: 3946

## Holdout performance

| model | horizon (min) | evaluated rows | MAE | RMSE |
| --- | ---: | ---: | ---: | ---: |
| Persistence (previous observation) | 10 | 3946 | 26.640 | 66.437 |
| Linear regression (weather + time + previous appliances) | 10 | 3946 | 26.211 | 61.354 |

## Reading of these numbers
- Against persistence, the regression lowers MAE from 26.640 to 26.211 and RMSE from 66.437 to 61.354.
- MAE improves by 0.429 Wh and RMSE by 5.083 Wh on the 3946 holdout steps.
- Robustness caveat: these are in-sample fitted, out-of-sample scored one-step-ahead forecasts from a single household, and the residuals in `results/md/residual_analysis.md` remain heavy-tailed.
