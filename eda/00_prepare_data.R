# WattWise -- EDA step 00: data preparation
# Downloads the UCI Appliances Energy Prediction dataset, validates its schema
# and timestamps, cleans conservatively, writes the chronological train/test
# split, and records the data-quality report.
#
# Run from the WattWise_Phase1 project root:
#   source("eda/00_prepare_data.R")
# or from a terminal:
#   Rscript eda/00_prepare_data.R

config <- list(
  source_file = "energydata_complete.csv",
  source_url = paste0(
    "https://archive.ics.uci.edu/static/public/374/",
    "appliances+energy+prediction.zip"
  ),
  train_fraction = 0.80,
  date_format = "%Y-%m-%d %H:%M:%S"
)

expected_columns <- c(
  "date", "Appliances", "lights",
  paste0("T", 1:9),
  paste0("RH_", 1:9),
  "T_out", "Press_mm_hg", "RH_out", "Windspeed",
  "Visibility", "Tdewpoint", "rv1", "rv2"
)

numeric_columns <- setdiff(expected_columns, "date")

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
      if (file.exists(file.path(path, "WattWise_Phase1.Rproj"))) {
        return(path)
      }
      if (dir.exists(file.path(path, "data")) &&
          dir.exists(file.path(path, "results"))) {
        return(path)
      }
      parent <- dirname(path)
      if (identical(parent, path)) break
      path <- parent
    }
  }
  
  stop(
    "Could not locate the WattWise project root from ", getwd(), ".",
    call. = FALSE
  )
}

download_uci_csv <- function(destination_dir) {
  dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
  
  zip_path <- tempfile(fileext = ".zip")
  on.exit(unlink(zip_path), add = TRUE)
  
  utils::download.file(
    config$source_url,
    destfile = zip_path,
    mode = "wb",
    method = "libcurl"
  )
  
  archive <- utils::unzip(zip_path, list = TRUE)
  matches <- archive$Name[
    basename(archive$Name) == config$source_file
  ]
  
  if (length(matches) != 1L) {
    stop(
      "Expected exactly one ", config$source_file,
      " in the UCI archive; found ", length(matches), "."
    )
  }
  
  utils::unzip(
    zip_path,
    files = matches,
    exdir = destination_dir,
    junkpaths = TRUE
  )
  
  destination <- file.path(destination_dir, config$source_file)
  if (!file.exists(destination) || file.info(destination)$size == 0) {
    stop("The UCI CSV was not downloaded or extracted correctly.")
  }
  
  destination
}

parse_timestamps <- function(x) {
  x <- trimws(as.character(x))
  parsed <- as.POSIXct(
    x,
    format = config$date_format,
    tz = "UTC"
  )
  
  invalid <- is.na(parsed)
  if (any(invalid)) {
    stop(
      "Invalid or missing timestamps at row(s): ",
      paste(head(which(invalid), 10L), collapse = ", "),
      ". Expected format: ", config$date_format
    )
  }
  
  # Require exact round-trip formatting so malformed dates are rejected.
  round_trip <- format(parsed, config$date_format, tz = "UTC")
  mismatch <- round_trip != x
  if (any(mismatch)) {
    stop(
      "Timestamps do not match the required format at row(s): ",
      paste(head(which(mismatch), 10L), collapse = ", ")
    )
  }
  
  parsed
}

write_timestamp_csv <- function(data, path) {
  output <- data
  output$date <- format(
    output$date,
    format = config$date_format,
    tz = "UTC"
  )
  utils::write.csv(output, path, row.names = FALSE, na = "")
}

run_pipeline <- function(force_download = FALSE) {
  project_root <- find_project_root()
  
  raw_dir <- file.path(project_root, "data", "raw")
  processed_dir <- file.path(project_root, "data", "processed")
  md_dir <- file.path(project_root, "results", "md")
  
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(md_dir, recursive = TRUE, showWarnings = FALSE)
  
  raw_path <- file.path(raw_dir, config$source_file)
  project_root_csv <- file.path(project_root, config$source_file)
  
  if (force_download && file.exists(raw_path)) {
    unlink(raw_path)
  }
  
  if (!file.exists(raw_path) && file.exists(project_root_csv)) {
    raw_path <- project_root_csv
    message("Using project-root dataset: ", raw_path)
  } else if (!file.exists(raw_path)) {
    raw_path <- download_uci_csv(raw_dir)
    message("Downloaded UCI dataset: ", raw_path)
  } else {
    message("Using local dataset: ", raw_path)
  }
  
  data <- utils::read.csv(
    raw_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("", "NA")
  )
  
  if (anyDuplicated(names(data))) {
    stop("Input contains duplicated column names.")
  }
  
  missing_columns <- setdiff(expected_columns, names(data))
  unexpected_columns <- setdiff(names(data), expected_columns)
  
  if (length(missing_columns) > 0L || length(unexpected_columns) > 0L) {
    stop(
      "Schema mismatch.\nMissing: ",
      paste(missing_columns, collapse = ", "),
      "\nUnexpected: ",
      paste(unexpected_columns, collapse = ", ")
    )
  }
  
  if (nrow(data) < 2L) {
    stop("Dataset must contain at least two rows.")
  }
  
  not_numeric <- numeric_columns[
    !vapply(data[numeric_columns], is.numeric, logical(1))
  ]
  if (length(not_numeric) > 0L) {
    stop("Expected numeric columns: ", paste(not_numeric, collapse = ", "))
  }
  
  parsed_date <- parse_timestamps(data$date)
  if (anyDuplicated(parsed_date)) {
    stop("Duplicate timestamps found; refusing ambiguous chronology.")
  }
  
  if (any(data$Appliances < 0, na.rm = TRUE)) {
    stop("Appliances contains negative target values.")
  }
  
  data$date <- parsed_date
  data <- data[order(data$date), , drop = FALSE]
  row.names(data) <- NULL
  
  numeric_matrix <- as.matrix(data[numeric_columns])
  finite_rows <- apply(numeric_matrix, 1L, function(x) all(is.finite(x)))
  keep_rows <- stats::complete.cases(data) & finite_rows
  
  cleaning <- list(
    raw_rows = nrow(data),
    excluded_rows = sum(!keep_rows),
    excluded_missing_rows = sum(!stats::complete.cases(data)),
    excluded_nonfinite_rows = sum(
      stats::complete.cases(data) & !finite_rows
    )
  )
  
  clean <- data[keep_rows, , drop = FALSE]
  row.names(clean) <- NULL
  
  if (nrow(clean) < 2L) {
    stop("Fewer than two usable rows remain after cleaning.")
  }
  
  n_train <- floor(nrow(clean) * config$train_fraction)
  n_train <- max(1L, min(n_train, nrow(clean) - 1L))
  
  train <- clean[seq_len(n_train), , drop = FALSE]
  test <- clean[seq.int(n_train + 1L, nrow(clean)), , drop = FALSE]
  
  if (!(max(train$date) < min(test$date))) {
    stop("Chronological split validation failed.")
  }
  
  write_timestamp_csv(
    clean,
    file.path(processed_dir, "energydata_clean.csv")
  )
  write_timestamp_csv(
    train,
    file.path(processed_dir, "train.csv")
  )
  write_timestamp_csv(
    test,
    file.path(processed_dir, "test.csv")
  )
  
  report <- c(
    "# Data quality report -- WattWise",
    "",
    paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste("Source file:", raw_path),
    "",
    "## Dimensions",
    paste("- Raw dimensions:", nrow(data), "rows x", ncol(data), "columns"),
    "- Schema exact match: TRUE",
    paste("- Raw missing cells:", sum(is.na(data))),
    paste("- Rows excluded by cleaning:", cleaning$excluded_rows),
    paste("  - rows with missing values:", cleaning$excluded_missing_rows),
    paste("  - rows with non-finite numeric values:",
          cleaning$excluded_nonfinite_rows),
    paste("- Clean rows:", nrow(clean)),
    "",
    "## Chronological split (80% train / 20% test, no shuffling)",
    paste("- Train rows:", nrow(train)),
    paste("- Train date range:",
          format(min(train$date), config$date_format, tz = "UTC"),
          "to",
          format(max(train$date), config$date_format, tz = "UTC")),
    paste("- Test rows:", nrow(test)),
    paste("- Test date range:",
          format(min(test$date), config$date_format, tz = "UTC"),
          "to",
          format(max(test$date), config$date_format, tz = "UTC")),
    paste("- Leakage boundary valid:", max(train$date) < min(test$date)),
    "",
    "## Cleaning policy",
    "Invalid or duplicate timestamps fail validation. Rows with missing or",
    "non-finite numeric values are excluded; no imputation or winsorising is",
    "applied. `rv1` and `rv2` are retained in the stored data but are never",
    "used as predictors."
  )
  
  writeLines(
    report,
    file.path(md_dir, "data_quality.md"),
    useBytes = TRUE
  )
  
  cat(paste(report, collapse = "\n"), "\n")
  invisible(list(clean = clean, train = train, test = test))
}

run_pipeline()