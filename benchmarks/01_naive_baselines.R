# WattWise -- benchmarks: naive forecast baselines
# Persistence (t-1) vs seasonal naive (same time previous day, 144 steps of
# 10 minutes). These are the reference numbers every linear model must beat.
#
# Run from the WattWise_Phase1 project root:
#   source("benchmarks/01_naive_baselines.R")
# or from a terminal:
#   Rscript benchmarks/01_naive_baselines.R

# Resolve the project root so the script runs from any working directory.
# Starts from the script's own location when run with Rscript, then from the
# working directory, and walks up until the project marker is found.
find_project_root <- function() {
  starts <- character(0)
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    starts <- c(starts, dirname(sub("^--file=", "", file_arg)))
  }
  starts <- c(starts, getwd())

  for (start in unique(starts)) {
    if (!dir.exists(start)) next
    path <- normalizePath(start, winslash = "/", mustWork = TRUE)
    repeat {
      if (file.exists(file.path(path, "WattWise_Phase1.Rproj"))) return(path)
      if (dir.exists(file.path(path, "data")) &&
          dir.exists(file.path(path, "results"))) return(path)
      parent <- dirname(path)
      if (identical(parent, path)) break
      path <- parent
    }
  }

  stop("Could not locate the WattWise project root from ", getwd(), ".",
       call. = FALSE)
}

project_root <- find_project_root()
source_path <- file.path(project_root, "data", "raw", "energydata_complete.csv")
csv_dir <- file.path(project_root, "results", "csv")
md_dir <- file.path(project_root, "results", "md")
season_lag <- 144L  # 24 hours at 10-minute intervals

if (!file.exists(source_path)) {
  stop("Missing source CSV: ", source_path)
}

data <- read.csv(
  source_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required <- c("date", "Appliances")
missing <- setdiff(required, names(data))
if (length(missing) > 0L) {
  stop("Missing columns: ", paste(missing, collapse = ", "))
}

# Parse timestamps with or without a time component.
raw_date <- trimws(as.character(data$date))
parsed_date <- as.POSIXct(
  raw_date,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"
)

date_only <- is.na(parsed_date) & grepl("^\\d{4}-\\d{2}-\\d{2}$", raw_date)
parsed_date[date_only] <- as.POSIXct(
  raw_date[date_only],
  format = "%Y-%m-%d",
  tz = "UTC"
)

if (anyNA(parsed_date)) {
  bad <- which(is.na(parsed_date))
  stop(
    "Unparsed timestamps at row(s): ",
    paste(head(bad, 10L), collapse = ", "),
    ". Check the source date format."
  )
}
if (anyDuplicated(parsed_date)) {
  stop("Duplicate timestamps found.")
}

data$date <- parsed_date
data <- data[order(data$date), , drop = FALSE]
row.names(data) <- NULL

if (!is.numeric(data$Appliances) ||
    anyNA(data$Appliances) ||
    any(!is.finite(data$Appliances)) ||
    any(data$Appliances < 0)) {
  stop("Appliances must be finite, non-missing, non-negative numeric values.")
}

# This comparison requires a complete 10-minute series.
interval_seconds <- as.numeric(diff(data$date), units = "secs")
if (any(interval_seconds != 600)) {
  bad_intervals <- which(interval_seconds != 600)
  stop(
    "The source series is not continuously sampled every 10 minutes. ",
    "Irregular interval count: ", length(bad_intervals),
    ". Do not impute date-only timestamps as midnight."
  )
}

# Use the same chronological 80/20 row split as Phase 1.
n_train <- floor(nrow(data) * 0.80)
test <- data[seq.int(n_train + 1L, nrow(data)), , drop = FALSE]
row.names(test) <- NULL

if (nrow(test) <= season_lag) {
  stop("Test partition must contain more than ", season_lag, " rows.")
}

# Compare both methods on identical test timestamps.
eval_rows <- seq.int(season_lag + 1L, nrow(test))
actual <- test$Appliances[eval_rows]

persistence_pred <- test$Appliances[eval_rows - 1L]
seasonal_pred <- test$Appliances[eval_rows - season_lag]

score <- function(model, prediction, observed, seasonal_lag_minutes = NA_integer_) {
  data.frame(
    model = model,
    horizon_minutes = 10L,
    seasonal_lag_minutes = seasonal_lag_minutes,
    evaluated_rows = length(observed),
    MAE = mean(abs(observed - prediction)),
    RMSE = sqrt(mean((observed - prediction)^2))
  )
}

metrics <- rbind(
  score(
    "Persistence (previous observation)",
    persistence_pred,
    actual
  ),
  score(
    "Seasonal naive (previous day)",
    seasonal_pred,
    actual,
    seasonal_lag_minutes = 1440L
  )
)

print(metrics)

dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(md_dir, showWarnings = FALSE, recursive = TRUE)
output_path <- file.path(csv_dir, "naive_baselines.csv")
write.csv(metrics, output_path, row.names = FALSE)

md_lines <- c(
  "# Naive forecast benchmarks -- WattWise",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Source script: `benchmarks/01_naive_baselines.R`",
  "",
  "One-step-ahead forecasts (10-minute horizon) scored on identical test",
  "timestamps from the chronological holdout period, so the two rows are",
  "directly comparable.",
  "",
  paste0(
    "| model | horizon (min) | seasonal lag (min) | evaluated rows | MAE | RMSE |"
  ),
  "| --- | ---: | ---: | ---: | ---: | ---: |",
  paste0(
    "| ", metrics$model,
    " | ", metrics$horizon_minutes,
    " | ", ifelse(is.na(metrics$seasonal_lag_minutes), "-", metrics$seasonal_lag_minutes),
    " | ", metrics$evaluated_rows,
    " | ", sprintf("%.3f", metrics$MAE),
    " | ", sprintf("%.3f", metrics$RMSE),
    " |"
  ),
  "",
  "## Reading of these numbers",
  paste0(
    "- Both baselines are scored on the same ", metrics$evaluated_rows[1L],
    " test timestamps."
  ),
  paste0(
    "- Persistence MAE ", sprintf("%.3f", metrics$MAE[1L]),
    " is lower than seasonal-naive MAE ", sprintf("%.3f", metrics$MAE[2L]),
    ", so copying the previous 10-minute reading is the stronger reference."
  ),
  paste0(
    "- Any regression model must improve on RMSE ",
    sprintf("%.3f", metrics$RMSE[1L]),
    " (persistence) to add value over the naive benchmark."
  )
)
writeLines(md_lines, file.path(md_dir, "naive_baselines.md"), useBytes = TRUE)

cat(
  "\nEvaluated timestamps:",
  format(min(test$date[eval_rows]), "%Y-%m-%d %H:%M:%S"),
  "to",
  format(max(test$date[eval_rows]), "%Y-%m-%d %H:%M:%S"),
  "\nSaved comparison to:", output_path,
  "\nSaved summary to:", file.path(md_dir, "naive_baselines.md"),
  "\n"
)