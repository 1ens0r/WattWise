# Naive forecast benchmarks -- WattWise

Generated: 2026-10-19 08:50:10 UTC
Source script: `benchmarks/01_naive_baselines.R`

One-step-ahead forecasts (10-minute horizon) scored on identical test
timestamps from the chronological holdout period, so the two rows are
directly comparable.

| model | horizon (min) | seasonal lag (min) | evaluated rows | MAE | RMSE |
| --- | ---: | ---: | ---: | ---: | ---: |
| Persistence (previous observation) | 10 | - | 3803 | 25.803 | 64.656 |
| Seasonal naive (previous day) | 10 | 1440 | 3803 | 51.420 | 108.013 |

## Reading of these numbers
- Both baselines are scored on the same 3803 test timestamps.
- Persistence MAE 25.803 is lower than seasonal-naive MAE 51.420, so copying the previous 10-minute reading is the stronger reference.
- Any regression model must improve on RMSE 64.656 (persistence) to add value over the naive benchmark.
