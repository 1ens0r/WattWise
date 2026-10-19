# Data quality report -- WattWise

Generated: 2026-10-19 08:50:10 UTC
Source file: C:/Users/sragh/Downloads/WattWise_Phase1/WattWise_Phase1/data/raw/energydata_complete.csv

## Dimensions
- Raw dimensions: 19735 rows x 29 columns
- Schema exact match: TRUE
- Raw missing cells: 0
- Rows excluded by cleaning: 0
  - rows with missing values: 0
  - rows with non-finite numeric values: 0
- Clean rows: 19735

## Chronological split (80% train / 20% test, no shuffling)
- Train rows: 15788
- Train date range: 2016-01-11 17:00:00 to 2016-04-30 08:10:00
- Test rows: 3947
- Test date range: 2016-04-30 08:20:00 to 2016-05-27 18:00:00
- Leakage boundary valid: TRUE

## Cleaning policy
Invalid or duplicate timestamps fail validation. Rows with missing or
non-finite numeric values are excluded; no imputation or winsorising is
applied. `rv1` and `rv2` are retained in the stored data but are never
used as predictors.
