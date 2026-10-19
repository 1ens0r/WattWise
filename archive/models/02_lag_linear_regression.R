# WattWise -- models: linear regression with six appliance lags
# Rolling one-step-ahead evaluation: at each test origin the six most recent
# observed appliance readings are available, together with weather and
# calendar time, to predict the next 10-minute reading.
#
# Run from the WattWise_Phase1 project root:
#   source("models/02_lag_linear_regression.R")
# or from a terminal:
#   Rscript models/02_lag_linear_regression.R

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
csv_dir <- file.path(project_root, "results", "csv")
md_dir <- file.path(project_root, "results", "md")

read_data <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  x$date <- as.POSIXct(
    x$date,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )
  if (anyNA(x$date)) stop("Could not parse all timestamps in ", path)
  x[order(x$date), , drop = FALSE]
}

all_data <- read_data(file.path(project_root, "data", "processed", "energydata_clean.csv"))
train <- read_data(file.path(project_root, "data", "processed", "train.csv"))
test <- read_data(file.path(project_root, "data", "processed", "test.csv"))

if (anyDuplicated(all_data$date)) stop("Duplicate timestamps found.")
if (max(train$date) >= min(test$date)) stop("Train and test periods overlap.")
if (any(as.numeric(diff(all_data$date), units = "secs") != 600)) {
  stop("Clean data is not regularly sampled every 10 minutes.")
}
if (!all(train$date %in% all_data$date) || !all(test$date %in% all_data$date)) {
  stop("Train/test files do not match the cleaned dataset.")
}

# Create one-step targets and lag features from the full ordered series.
# An example at row i uses Appliances[i] ... Appliances[i-5] to predict
# Appliances[i+1]. Keep training targets in train and test targets in test.
n <- nrow(all_data)
if (n < 8L) stop("Not enough rows to create six lag features.")

origin_idx <- 6L:(n - 1L)
target_idx <- origin_idx + 1L

make_features <- function(i) {
  d <- all_data
  hour <- as.integer(format(d$date[i], "%H", tz = "UTC")) +
    as.integer(format(d$date[i], "%M", tz = "UTC")) / 60
  weekday <- as.integer(format(d$date[i], "%u", tz = "UTC"))
  
  out <- data.frame(
    appliance_lag1 = d$Appliances[i],
    appliance_lag2 = d$Appliances[i - 1L],
    appliance_lag3 = d$Appliances[i - 2L],
    appliance_lag4 = d$Appliances[i - 3L],
    appliance_lag5 = d$Appliances[i - 4L],
    appliance_lag6 = d$Appliances[i - 5L],
    T_out = d$T_out[i],
    RH_out = d$RH_out[i],
    Windspeed = d$Windspeed[i],
    Visibility = d$Visibility[i],
    Press_mm_hg = d$Press_mm_hg[i],
    hour_sin = sin(2 * pi * hour / 24),
    hour_cos = cos(2 * pi * hour / 24),
    weekday_sin = sin(2 * pi * weekday / 7),
    weekday_cos = cos(2 * pi * weekday / 7)
  )
  out
}

features <- do.call(rbind, lapply(origin_idx, make_features))
target_date <- all_data$date[target_idx]
target_value <- all_data$Appliances[target_idx]

train_rows <- target_date <= max(train$date)
test_rows <- target_date >= min(test$date) & target_date <= max(test$date)

if (!any(train_rows) || !any(test_rows)) stop("Could not form train/test examples.")

train_x <- features[train_rows, , drop = FALSE]
test_x <- features[test_rows, , drop = FALSE]
train_y <- target_value[train_rows]
test_y <- target_value[test_rows]
test_target_dates <- target_date[test_rows]

# Verify all model training targets precede the test period.
if (max(target_date[train_rows]) >= min(test_target_dates)) {
  stop("Training targets overlap the test period.")
}

fit_data <- train_x
fit_data$target <- train_y
fit <- lm(target ~ ., data = fit_data)

regression_prediction <- as.numeric(predict(fit, newdata = test_x))
persistence_prediction <- test_x$appliance_lag1

score <- function(model, prediction, observed) {
  data.frame(
    model = model,
    horizon_minutes = 10L,
    evaluated_rows = length(observed),
    MAE = mean(abs(observed - prediction)),
    RMSE = sqrt(mean((observed - prediction)^2))
  )
}

metrics <- rbind(
  score("Persistence (previous observation)",
        persistence_prediction, test_y),
  score("Linear regression (six appliance lags + weather + time)",
        regression_prediction, test_y)
)

print(metrics)

errors <- data.frame(
  date = format(test_target_dates, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
  actual = test_y,
  persistence_prediction = persistence_prediction,
  regression_prediction = regression_prediction,
  persistence_absolute_error = abs(test_y - persistence_prediction),
  regression_absolute_error = abs(test_y - regression_prediction)
)

dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(md_dir, showWarnings = FALSE, recursive = TRUE)
metrics_path <- file.path(csv_dir, "lag_regression_metrics.csv")
errors_path <- file.path(csv_dir, "lag_regression_errors.csv")
write.csv(metrics, metrics_path, row.names = FALSE)
write.csv(errors, errors_path, row.names = FALSE)

persistence_mae <- metrics$MAE[1L]
persistence_rmse <- metrics$RMSE[1L]
model_mae <- metrics$MAE[2L]
model_rmse <- metrics$RMSE[2L]

md_lines <- c(
  "# Linear regression with six appliance lags -- WattWise",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Source script: `models/02_lag_linear_regression.R`",
  "",
  "## Setup",
  "- Target: `Appliances` at time *t+1* (watt-hours), one 10-minute step ahead.",
  paste("- Predictors at time *t*: six appliance lags, outdoor weather",
        "(`T_out`, `RH_out`, `Windspeed`, `Visibility`, `Press_mm_hg`) and",
        "sine/cosine terms for hour of day and weekday."),
  "- Estimator: ordinary least squares via `lm()` fitted on training targets only.",
  paste("- Evaluation: rolling one-step-ahead on the chronological holdout,",
        metrics$evaluated_rows[2L], "steps."),
  "",
  "## Holdout performance",
  "",
  "| model | horizon (min) | evaluated rows | MAE | RMSE |",
  "| --- | ---: | ---: | ---: | ---: |",
  paste0(
    "| ", metrics$model,
    " | ", metrics$horizon_minutes,
    " | ", metrics$evaluated_rows,
    " | ", sprintf("%.3f", metrics$MAE),
    " | ", sprintf("%.3f", metrics$RMSE),
    " |"
  ),
  "",
  "## Reading of these numbers",
  paste0(
    "- Against persistence at the same origins, MAE moves from ",
    sprintf("%.3f", persistence_mae), " to ", sprintf("%.3f", model_mae),
    " and RMSE from ", sprintf("%.3f", persistence_rmse), " to ",
    sprintf("%.3f", model_rmse), "."
  ),
  paste0(
    "- Per-row actual, prediction and absolute error for both models are in",
    " `results/csv/lag_regression_errors.csv` (", nrow(errors), " rows)."
  ),
  paste0(
    "- RMSE remains far above MAE for both models, which indicates a",
    " small number of large spikes dominates the squared-error loss."
  )
)
writeLines(md_lines, file.path(md_dir, "lag_regression.md"), useBytes = TRUE)

cat(
  "\nSaved metrics:", metrics_path,
  "\nSaved row errors:", errors_path,
  "\nSaved summary:", file.path(md_dir, "lag_regression.md"), "\n"
)