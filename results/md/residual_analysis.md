# Residual and error analysis -- WattWise

Generated: 2026-10-19 08:10:27 UTC
Source script: `models/03_residual_analysis.R`

## Setup
- Model: the linear regression of `models/01_linear_regression.R`
  (previous appliances, outdoor weather, hour and weekday).
- Scope: 3946 scored holdout steps. Absolute errors are Actual - Prediction, in watt-hours.

## Error by hour of day

| hour | rows | persistence MAE | persistence RMSE | regression MAE | regression RMSE |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 162 | 7.284 | 10.244 | 6.966 | 9.074 |
| 1 | 162 | 6.605 | 9.196 | 7.004 | 8.658 |
| 2 | 162 | 8.951 | 11.358 | 8.905 | 10.625 |
| 3 | 162 | 7.963 | 10.801 | 8.334 | 10.218 |
| 4 | 162 | 6.235 | 8.784 | 6.574 | 8.285 |
| 5 | 162 | 15.926 | 38.329 | 14.635 | 34.117 |
| 6 | 162 | 20.864 | 47.635 | 18.629 | 42.109 |
| 7 | 162 | 38.210 | 91.054 | 35.121 | 85.219 |
| 8 | 165 | 56.788 | 105.207 | 50.975 | 93.142 |
| 9 | 168 | 45.000 | 98.410 | 43.675 | 97.501 |
| 10 | 168 | 43.929 | 91.716 | 38.459 | 79.146 |
| 11 | 168 | 37.798 | 94.878 | 39.446 | 89.462 |
| 12 | 168 | 31.726 | 71.034 | 33.433 | 66.381 |
| 13 | 168 | 30.714 | 74.434 | 34.374 | 66.691 |
| 14 | 168 | 21.845 | 59.577 | 26.990 | 57.262 |
| 15 | 168 | 18.095 | 39.940 | 23.600 | 38.245 |
| 16 | 168 | 39.464 | 83.456 | 42.628 | 81.441 |
| 17 | 168 | 69.167 | 122.372 | 64.948 | 119.521 |
| 18 | 163 | 53.681 | 102.953 | 46.643 | 83.145 |
| 19 | 162 | 15.679 | 26.364 | 16.313 | 24.288 |
| 20 | 162 | 23.272 | 54.337 | 21.306 | 49.622 |
| 21 | 162 | 14.691 | 32.069 | 14.750 | 28.179 |
| 22 | 162 | 10.000 | 12.910 | 9.884 | 12.211 |
| 23 | 162 | 11.111 | 31.972 | 10.739 | 28.549 |

## Ten largest regression misses

| date (UTC) | hour | actual | persistence pred | regression pred | persistence abs err | regression abs err |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 2016-05-06 07:30:00 | 7 | 670.0 | 70.0 | 69.0 | 600.0 | 601.0 |
| 2016-05-27 09:40:00 | 9 | 580.0 | 50.0 | 61.4 | 530.0 | 518.6 |
| 2016-05-16 17:40:00 | 17 | 630.0 | 110.0 | 113.9 | 520.0 | 516.1 |
| 2016-05-26 11:40:00 | 11 | 600.0 | 100.0 | 102.0 | 500.0 | 498.0 |
| 2016-05-13 11:10:00 | 11 | 600.0 | 110.0 | 110.5 | 490.0 | 489.5 |
| 2016-05-13 07:50:00 | 7 | 560.0 | 80.0 | 78.6 | 480.0 | 481.4 |
| 2016-05-21 09:40:00 | 9 | 790.0 | 390.0 | 310.1 | 400.0 | 479.9 |
| 2016-05-21 11:10:00 | 11 | 590.0 | 110.0 | 111.5 | 480.0 | 478.5 |
| 2016-04-30 12:20:00 | 12 | 570.0 | 80.0 | 94.1 | 490.0 | 475.9 |
| 2016-05-07 18:10:00 | 18 | 560.0 | 80.0 | 89.9 | 480.0 | 470.1 |

## Reading of these numbers
- Error is clearly time-of-day dependent: regression MAE peaks in hour 17 (64.948 Wh) and is lowest in hour 4 (6.574 Wh).
- The largest single regression error is 601.0 Wh on 2016-05-06 07:30:00, a spike both models miss; persistence and the regression fail on the
  same events, which is the signature of unobserved appliance switching rather than a modelling artefact.
- Errors are not independent over time: neighbouring 10-minute steps share the same spikes, so the usual i.i.d. inference caveat applies.
