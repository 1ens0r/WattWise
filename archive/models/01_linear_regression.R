# WattWise -- models: linear regression (weather, calendar time, lag-1)
# Ordinary least squares one-step-ahead prediction of Appliances at t+1 from
# the observation at t: previous appliances, outdoor weather, hour and weekday.
# Benchmarked against persistence on the chronological holdout.
#
# Run from the WattWise_Phase1 project root:
#   source("models/01_linear_regression.R")
# or from a terminal:
#   Rscript models/01_linear_regression.R

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
clean_path <- file.path(project_root, "data", "processed", "energydata_clean.csv")
train_path <- file.path(project_root, "data", "processed", "train.csv")
test_path  <- file.path(project_root, "data", "processed", "test.csv")
csv_dir <- file.path(project_root, "results", "csv")
md_dir <- file.path(project_root, "results", "md")

for (path in c(clean_path, train_path, test_path)) {
  if (!file.exists(path)) stop("Missing file: ", path, ". Run Phase 1 first.")
}

read_wattwise <- function(path) {
  x <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  x$date <- as.POSIXct(
    x$date,
    format = "%Y-%m-%d %H:%M:%S",
    tz = "UTC"
  )
  if (anyNA(x$date)) stop("Could not parse all timestamps in: ", path)
  x[order(x$date), , drop = FALSE]
}

data <- read_wattwise(clean_path)
train <- read_wattwise(train_path)
test <- read_wattwise(test_path)

required <- c(
  "date", "Appliances", "T_out", "RH_out",
  "Windspeed", "Visibility", "Press_mm_hg"
)
missing <- setdiff(required, names(data))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

if (anyDuplicated(data$date) || anyDuplicated(train$date) || anyDuplicated(test$date)) {
  stop("Duplicate timestamps found.")
}

# Ensure Phase 1 partitions correspond to the cleaned full series.
if (max(train$date) >= min(test$date)) stop("Train/test periods overlap.")
if (!all(train$date %in% data$date) || !all(test$date %in% data$date)) {
  stop("Train/test files do not match energydata_clean.csv.")
}

# Build features at time t and target Appliances at t+1.
# Use only pairs whose two observations are within the same partition.
make_pairs <- function(partition) {
  if (nrow(partition) < 2L) stop("Partition needs at least two rows.")
  dt <- as.numeric(diff(partition$date), units = "secs")
  if (any(dt != 600)) stop("Partition is not uniformly sampled every 10 minutes.")
  
  origin <- partition[-nrow(partition), , drop = FALSE]
  target <- partition$Appliances[-1L]
  
  hour <- as.integer(format(origin$date, "%H", tz = "UTC")) +
    as.integer(format(origin$date, "%M", tz = "UTC")) / 60
  
  predictors <- data.frame(
    previous_appliances = origin$Appliances,
    T_out = origin$T_out,
    RH_out = origin$RH_out,
    Windspeed = origin$Windspeed,
    Visibility = origin$Visibility,
    Press_mm_hg = origin$Press_mm_hg,
    hour_sin = sin(2 * pi * hour / 24),
    hour_cos = cos(2 * pi * hour / 24),
    weekday_sin = sin(2 * pi * as.integer(format(origin$date, "%u", tz = "UTC")) / 7),
    weekday_cos = cos(2 * pi * as.integer(format(origin$date, "%u", tz = "UTC")) / 7)
  )
  
  list(origin_date = origin$date, x = predictors, y = target)
}

train_pairs <- make_pairs(train)
test_pairs <- make_pairs(test)

if (anyNA(train_pairs$x) || anyNA(test_pairs$x) ||
    anyNA(train_pairs$y) || anyNA(test_pairs$y)) {
  stop("Missing predictor or target values found.")
}

model_data <- train_pairs$x
model_data$target_next_appliances <- train_pairs$y

fit <- lm(target_next_appliances ~ ., data = model_data)

regression_prediction <- as.numeric(
  predict(fit, newdata = test_pairs$x)
)
persistence_prediction <- test_pairs$x$previous_appliances
actual <- test_pairs$y

score <- function(name, prediction, observed) {
  data.frame(
    model = name,
    horizon_minutes = 10L,
    evaluated_rows = length(observed),
    MAE = mean(abs(observed - prediction)),
    RMSE = sqrt(mean((observed - prediction)^2))
  )
}

metrics <- rbind(
  score("Persistence (previous observation)", persistence_prediction, actual),
  score("Linear regression (weather + time + previous appliances)",
        regression_prediction, actual)
)

print(metrics)

test_start <- format(min(test_pairs$origin_date + 600), "%Y-%m-%d %H:%M:%S")
test_end <- format(max(test_pairs$origin_date + 600), "%Y-%m-%d %H:%M:%S")
cat("\nTest target dates:", test_start, "to", test_end, "\n")

dir.create(csv_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(md_dir, showWarnings = FALSE, recursive = TRUE)
output_path <- file.path(csv_dir, "linear_regression_vs_persistence.csv")
write.csv(metrics, output_path, row.names = FALSE)

persistence_mae <- metrics$MAE[1L]
persistence_rmse <- metrics$RMSE[1L]
model_mae <- metrics$MAE[2L]
model_rmse <- metrics$RMSE[2L]

md_lines <- c(
  "# Linear regression (weather + calendar time + lag-1) -- WattWise",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Source script: `models/01_linear_regression.R`",
  "",
  "## Setup",
  paste("- Target: `Appliances` at time *t+1* (watt-hours)."),
  paste("- Predictors: measurements at time *t*:",
        paste0("`", paste(names(test_pairs$x), collapse = "`, `"), "`.")),
  paste("- Estimator: ordinary least squares via `lm()` on the training",
        "partition only; the test partition is never used for fitting."),
  paste("- Evaluation period:", test_start, "to", test_end),
  paste("- Evaluated rows:", metrics$evaluated_rows[2L]),
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
    "- Against persistence, the regression lowers MAE from ",
    sprintf("%.3f", persistence_mae), " to ", sprintf("%.3f", model_mae),
    " and RMSE from ", sprintf("%.3f", persistence_rmse), " to ",
    sprintf("%.3f", model_rmse), "."
  ),
  paste0(
    "- MAE improves by ", sprintf("%.3f", persistence_mae - model_mae),
    " Wh and RMSE by ", sprintf("%.3f", persistence_rmse - model_rmse),
    " Wh on the ", metrics$evaluated_rows[2L], " holdout steps."
  ),
  paste0(
    "- Robustness caveat: these are in-sample fitted, out-of-sample scored",
    " one-step-ahead forecasts from a single household, and the residuals in",
    " `results/md/residual_analysis.md` remain heavy-tailed."
  )
)
writeLines(md_lines, file.path(md_dir, "linear_regression.md"), useBytes = TRUE)

cat("Saved metrics to:", output_path, "\n")
cat("Saved summary to:", file.path(md_dir, "linear_regression.md"), "\n")