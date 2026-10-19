# Linear regression with six appliance lags -- WattWise

Generated: 2026-10-19 08:10:24 UTC
Source script: `models/02_lag_linear_regression.R`

## Setup
- Target: `Appliances` at time *t+1* (watt-hours), one 10-minute step ahead.
- Predictors at time *t*: six appliance lags, outdoor weather (`T_out`, `RH_out`, `Windspeed`, `Visibility`, `Press_mm_hg`) and sine/cosine terms for hour of day and weekday.
- Estimator: ordinary least squares via `lm()` fitted on training targets only.
- Evaluation: rolling one-step-ahead on the chronological holdout, 3947 steps.

## Holdout performance

| model | horizon (min) | evaluated rows | MAE | RMSE |
| --- | ---: | ---: | ---: | ---: |
| Persistence (previous observation) | 10 | 3947 | 26.651 | 66.438 |
| Linear regression (six appliance lags + weather + time) | 10 | 3947 | 26.365 | 60.114 |

## Reading of these numbers
- Against persistence at the same origins, MAE moves from 26.651 to 26.365 and RMSE from 66.438 to 60.114.
- Per-row actual, prediction and absolute error for both models are in `results/csv/lag_regression_errors.csv` (3947 rows).
- RMSE remains far above MAE for both models, which indicates a small number of large spikes dominates the squared-error loss.
