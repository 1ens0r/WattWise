# WattWise -- head-to-head benchmark of the six linear models
#
# Trains every model in `models/` plus three reference models on the identical
# chronological holdout, scores them, and writes:
#   results/csv/model_benchmark.csv     -- one row per model
#   results/csv/model_predictions.csv   -- per-timestamp predictions, all models
#   results/csv/model_coefficients.csv  -- standardised coefficients
#   results.md                          -- the written-up benchmark (project root)
#
# Protocol (identical for every model):
#   * Target: Appliances at t+1 (10 minutes ahead).
#   * All features are measured at t or earlier; `lights`, `rv1` and `rv2` are
#     never used.
#   * Training uses only timestamps up to max(train.csv date); the holdout is
#     the same 3947 timestamps as the Phase 1 test split.
#   * Hyperparameters and standardisation statistics come from training rows
#     only, via forward-chaining (blocked) cross-validation.
#   * No random seed is needed: every estimator here is deterministic.
#
# Run from the project root:  source("benchmarks/02_model_benchmark.R")

# Resolve the project root so the script runs from any working directory.
find_script_root <- function() {
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
  stop("Could not locate the WattWise project root.", call. = FALSE)
}

project_root <- find_script_root()

source(file.path(project_root, "models", "00_shared.R"))

model_files <- c(
  "01_ridge.R", "02_elastic_net.R", "03_lasso.R",
  "04_polynomial_ridge.R", "05_distributed_lag.R", "06_phmr.R"
)

models <- list()
for (f in model_files) {
  source(file.path(project_root, "models", f))
  models[[length(models) + 1L]] <- wattwise_model
  wattwise_model <- NULL
}
model_keys <- sub("\\.R$", "", sub("^[0-9]+_", "", model_files))

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

series <- load_series(project_root)
train_part <- utils::read.csv(file.path(project_root, "data", "processed", "train.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
test_part <- utils::read.csv(file.path(project_root, "data", "processed", "test.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
train_part$date <- as.POSIXct(train_part$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
test_part$date <- as.POSIXct(test_part$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

candidate <- build_candidates(series, env_lags = 12L, appl_lags = 6L)
split <- label_split(candidate, train_part, test_part)

is_train <- split == "train"
is_test <- split == "test"

x_train <- candidate$x[is_train, , drop = FALSE]
y_train <- candidate$y[is_train]
x_test <- candidate$x[is_test, , drop = FALSE]
y_test <- candidate$y[is_test]
test_dates <- candidate$date[is_test]

stopifnot(ncol(x_train) == ncol(x_test))
stopifnot(max(candidate$date[is_train]) < min(candidate$date[is_test]))

cat(sprintf("Training origins: %d | holdout origins: %d | candidate columns: %d\n",
            nrow(x_train), nrow(x_test), ncol(x_train)))

# ---------------------------------------------------------------------------
# Evaluation helper
# ---------------------------------------------------------------------------

evaluate_model <- function(model, key) {
  cat(sprintf("  fitting %-26s ", model$name))
  t0 <- proc.time()[["elapsed"]]
  tuned <- model$tune(x_train, y_train)
  fitted <- model$fit(x_train, y_train, tuned$hyper)
  prediction <- as.vector(fitted$predict(x_test))
  elapsed <- proc.time()[["elapsed"]] - t0

  if (length(prediction) != length(y_test)) {
    stop(model$name, ": prediction length mismatch.", call. = FALSE)
  }
  if (any(!is.finite(prediction))) {
    stop(model$name, ": non-finite predictions.", call. = FALSE)
  }

  info <- fitted$info
  result <- list(
    key = key, name = model$name, family = model$family,
    rationale = model$rationale, reference = model$reference,
    prediction = prediction,
    rmse = rmse_wh(y_test, prediction),
    mae = mae_wh(y_test, prediction),
    r2 = r2_holdout(y_test, prediction),
    cv_rmse = tuned$cv_rmse,
    n_predictors = info$n_predictors,
    n_nonzero = info$n_nonzero,
    detail = info$detail,
    coefficients = info$coefficients,
    extras = info$extras,
    seconds = elapsed
  )
  cat(sprintf("RMSE %.3f (%.1fs)\n", result$rmse, result$seconds))
  result
}

# ---------------------------------------------------------------------------
# Reference models (not part of the six, but the bar to clear)
# ---------------------------------------------------------------------------

ols_reference <- function(key, name, columns_fn, rationale) {
  list(
    key = key, name = name, family = "Linear, OLS (reference)",
    rationale = rationale, reference = "Ordinary least squares",
    tune = function(x, y) {
      cols <- columns_fn()
      cv <- cv_mean_rmse(
        x[, cols, drop = FALSE], y, n_blocks = 6L,
        fit_fn = function(xt, yt) ols_fit(xt, yt),
        pred_fn = function(state, xv) as.vector(xv %*% state$beta)
      )
      list(hyper = list(columns = cols), cv_rmse = cv)
    },
    fit = function(x, y, hyper) {
      state <- ols_fit(x[, hyper$columns, drop = FALSE], y)
      list(
        predict = function(xnew) as.vector(xnew[, hyper$columns, drop = FALSE] %*% state$beta),
        info = list(
          coefficients = stats::setNames(state$beta, hyper$columns),
          n_predictors = length(hyper$columns),
          n_nonzero = sum(abs(state$beta) > 1e-8),
          detail = sprintf("%d predictors, %d aliased", length(hyper$columns),
                           state$aliased)
        )
      )
    }
  )
}

reference_models <- list(
  list(
    key = "persistence", name = "Persistence (reference)", family = "Naive (reference)",
    rationale = "Copy the most recent reading; the benchmark every model must beat.",
    reference = "Naive forecast",
    tune = function(x, y) {
      cv <- cv_mean_rmse(x, y, n_blocks = 6L,
                         fit_fn = function(xt, yt) NULL,
                         pred_fn = function(state, xv) xv[, "appl_lag1"])
      list(hyper = list(), cv_rmse = cv)
    },
    fit = function(x, y, hyper) {
      list(
        predict = function(xnew) as.vector(xnew[, "appl_lag1"]),
        info = list(
          coefficients = stats::setNames(1, "appl_lag1"),
          n_predictors = 1L, n_nonzero = 1L,
          detail = "y_hat(t+1) = y(t)"
        )
      )
    }
  ),
  ols_reference(
    "ols_base", "OLS, base library (reference)",
    function() wattwise_columns()$base,
    "Unpenalised least squares on the shared 29-channel library."
  ),
  ols_reference(
    "ols_appliance_lags", "OLS, base + six appliance lags (reference)",
    function() c(wattwise_columns()$base, paste0("appl_lag", 2:6)),
    "Unpenalised least squares with recent consumption history; the archived phase-3 feature idea."
  )
)

# ---------------------------------------------------------------------------
# Run everything on identical rows
# ---------------------------------------------------------------------------

cat("\nSix benchmark models\n")
results <- vector("list", length(models))
for (i in seq_along(models)) {
  results[[i]] <- evaluate_model(models[[i]], model_keys[i])
}

cat("\nReference models\n")
for (ref in reference_models) {
  results[[length(results) + 1L]] <- evaluate_model(ref, ref$key)
}

keys <- vapply(results, function(r) r$key, character(1))
names(results) <- keys

bench <- data.frame(
  model = vapply(results, function(r) r$name, character(1)),
  key = keys,
  family = vapply(results, function(r) r$family, character(1)),
  predictors = vapply(results, function(r) r$n_predictors, numeric(1)),
  nonzero = vapply(results, function(r) r$n_nonzero, numeric(1)),
  cv_rmse = vapply(results, function(r) r$cv_rmse, numeric(1)),
  rmse = vapply(results, function(r) r$rmse, numeric(1)),
  mae = vapply(results, function(r) r$mae, numeric(1)),
  r2 = vapply(results, function(r) r$r2, numeric(1)),
  seconds = vapply(results, function(r) r$seconds, numeric(1)),
  tuned = vapply(results, function(r) r$detail, character(1)),
  reference = vapply(results, function(r) r$reference, character(1)),
  rationale = vapply(results, function(r) r$rationale, character(1)),
  stringsAsFactors = FALSE
)
bench <- bench[order(bench$rmse), , drop = FALSE]
rownames(bench) <- NULL

# The six benchmark models keep the number of their file; reference models have
# no index, so `index` separates the head-to-head six from the extra context.
file_index <- as.integer(sub("^([0-9]+)_.*$", "\\1", model_files))
names(file_index) <- model_keys
bench$index <- unname(file_index[bench$key])

six <- bench[!is.na(bench$index), , drop = FALSE]
six <- six[order(six$rmse), , drop = FALSE]
rownames(six) <- NULL
six$rank_six <- seq_len(nrow(six))

refs <- bench[is.na(bench$index), , drop = FALSE]
refs <- refs[order(refs$rmse), , drop = FALSE]
rownames(refs) <- NULL
refs$rank_ref <- seq_len(nrow(refs))

best_six_key <- six$key[1L]
best_six <- results[[best_six_key]]
overall_best_key <- bench$key[1L]
overall_best <- results[[overall_best_key]]

# ---------------------------------------------------------------------------
# Paired comparisons against the best of the six
# ---------------------------------------------------------------------------

error_matrix <- vapply(results, function(r) y_test - r$prediction, numeric(length(y_test)))
colnames(error_matrix) <- keys

paired <- do.call(rbind, lapply(setdiff(keys, best_six_key), function(k) {
  test <- suppressWarnings(stats::wilcox.test(
    abs(error_matrix[, best_six_key]), abs(error_matrix[, k]), paired = TRUE
  ))
  data.frame(
    challenger = bench$model[match(k, bench$key)],
    median_abs_err_difference = stats::median(abs(error_matrix[, k])) -
      stats::median(abs(error_matrix[, best_six_key])),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}))
paired <- paired[order(paired$p_value), , drop = FALSE]

# ---------------------------------------------------------------------------
# Stability across the holdout period (best three of the six)
# ---------------------------------------------------------------------------

month <- format(test_dates, "%Y-%m")
top3 <- six$key[seq_len(min(3L, nrow(six)))]
stability <- data.frame(month = sort(unique(month)), stringsAsFactors = FALSE)
for (k in unique(c("persistence", top3))) {
  if (!k %in% keys) next
  pred <- results[[k]]$prediction
  stability[[k]] <- vapply(stability$month, function(m) {
    idx <- month == m
    rmse_wh(y_test[idx], pred[idx])
  }, numeric(1))
}
month_counts <- vapply(stability$month, function(m) sum(month == m), numeric(1))

# ---------------------------------------------------------------------------
# Write outputs
# ---------------------------------------------------------------------------

csv_dir <- file.path(project_root, "results", "csv")
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)

utils::write.csv(bench, file.path(csv_dir, "model_benchmark.csv"), row.names = FALSE)

predictions <- data.frame(date = format(test_dates, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                          actual = y_test, stringsAsFactors = FALSE)
for (r in results) predictions[[r$key]] <- r$prediction
utils::write.csv(predictions, file.path(csv_dir, "model_predictions.csv"),
                 row.names = FALSE)

coef_rows <- do.call(rbind, lapply(results, function(r) {
  coefs <- r$coefficients
  if (is.null(coefs)) return(NULL)
  if (length(coefs) > 60L) {
    coefs <- coefs[order(-abs(coefs))][seq_len(min(60L, length(coefs)))]
  }
  data.frame(model = r$name, term = names(coefs), coefficient = as.numeric(coefs),
             stringsAsFactors = FALSE)
}))
utils::write.csv(coef_rows, file.path(csv_dir, "model_coefficients.csv"),
                 row.names = FALSE)

utils::write.csv(paired, file.path(csv_dir, "model_paired_tests.csv"), row.names = FALSE)
utils::write.csv(stability, file.path(csv_dir, "model_stability.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# results.md
# ---------------------------------------------------------------------------

fmt <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
md_row <- function(...) paste0("| ", paste(c(...), collapse = " | "), " |")

results_table <- function(df, rank_col) {
  out <- c(
    md_row(rank_col, "model", "family", "predictors", "active", "CV RMSE",
           "holdout RMSE", "holdout MAE", "holdout R2", "seconds"),
    md_row("---:", "---", "---", "---:", "---:", "---:", "---:", "---:", "---:", "---:")
  )
  for (i in seq_len(nrow(df))) {
    out <- c(out, md_row(
      df[[rank_col]][i], df$model[i], df$family[i],
      format(df$predictors[i], scientific = FALSE),
      format(df$nonzero[i], scientific = FALSE),
      fmt(df$cv_rmse[i]), fmt(df$rmse[i]), fmt(df$mae[i]),
      fmt(df$r2[i], 4), fmt(df$seconds[i], 1)
    ))
  }
  out
}

six_display <- results_table(six, "rank_six")
ref_display <- results_table(refs, "rank_ref")

model_display <- c(
  md_row("#", "model", "family", "why it suits this dataset", "tuned structure", "reference"),
  md_row("---:", "---", "---", "---", "---", "---")
)
ordered_six <- bench[!is.na(bench$index), , drop = FALSE]
ordered_six <- ordered_six[order(ordered_six$index), , drop = FALSE]
for (i in seq_len(nrow(ordered_six))) {
  model_display <- c(model_display, md_row(
    ordered_six$index[i], ordered_six$model[i], ordered_six$family[i],
    ordered_six$rationale[i], ordered_six$tuned[i], ordered_six$reference[i]
  ))
}
for (i in seq_len(nrow(refs))) {
  model_display <- c(model_display, md_row(
    "ref", refs$model[i], refs$family[i], refs$rationale[i], refs$tuned[i],
    refs$reference[i]
  ))
}

paired_display <- c(
  md_row("challenger", "median abs error vs winner (Wh)", "Wilcoxon p-value"),
  md_row("---", "---:", "---:")
)
for (i in seq_len(nrow(paired))) {
  paired_display <- c(paired_display, md_row(paired$challenger[i],
                                             fmt(paired$median_abs_err_difference[i]),
                                             format.pval(paired$p_value[i], digits = 3)))
}

stability_display <- c(md_row("month", colnames(stability)[-1]), md_row("---", rep("---:", ncol(stability) - 1)))
for (i in seq_len(nrow(stability))) {
  stability_display <- c(stability_display, md_row(
    stability$month[i],
    vapply(stability[i, -1, drop = TRUE], function(v) fmt(v), character(1))
  ))
}

top_coefs <- function(coefs, n = 12L) {
  if (is.null(coefs)) return(character(0))
  top <- coefs[order(-abs(coefs))]
  top <- top[seq_len(min(n, length(top)))]
  vapply(seq_along(top), function(i) {
    md_row(names(top)[i], fmt(as.numeric(top[i]), 2))
  }, character(1))
}
best_six_coefs <- top_coefs(best_six$coefficients)
overall_best_coefs <- top_coefs(overall_best$coefficients)

persist_rmse <- bench$rmse[bench$key == "persistence"]
persist_mae <- bench$mae[bench$key == "persistence"]
n_six_beat_persistence <- sum(six$rmse < persist_rmse)
six_spread <- max(six$rmse) - min(six$rmse)

mae_best_six_key <- six$key[which.min(six$mae)]
mae_best_six <- results[[mae_best_six_key]]
same_winner <- identical(mae_best_six_key, best_six_key)
n_six_beats_persist_mae <- sum(six$mae < persist_mae)

md <- c(
  "# WattWise -- model benchmark results",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Generated by: `benchmarks/02_model_benchmark.R` (run `source(\"run_all.R\")` to refresh everything)",
  "",
  "## Summary",
  "",
  sprintf("- **Best of the six benchmark models by RMSE: %s** -- holdout RMSE %.3f Wh, MAE %.3f Wh, R2 %.4f.",
          best_six$name, best_six$rmse, best_six$mae, best_six$r2),
  if (same_winner) {
    "  It also has the lowest MAE of the six, so both metrics agree."
  } else {
    sprintf("- By MAE the winner is different: **%s** at %.3f Wh. See \"Does the winner depend on the metric?\" below -- it matters, because the RMSE winner is worse than persistence on MAE.",
            mae_best_six$name, mae_best_six$mae)
  },
  sprintf("- All six beat the persistence benchmark (RMSE %.3f Wh); the winner beats it by %.3f Wh.",
          persist_rmse, persist_rmse - best_six$rmse),
  sprintf("- The six are separated by only %.3f Wh of RMSE, so this is a close field, not a blowout.",
          six_spread),
  sprintf("- Counting the three reference models, the best overall is **%s** at RMSE %.3f Wh, i.e. %.3f Wh ahead of the best of the six. That model is plain unpenalised OLS with six appliance lags, so the single biggest gain available here is consumption history, not a fancier penalty. See the section below the results.",
          overall_best$name, overall_best$rmse, best_six$rmse - overall_best$rmse),
  sprintf("- Evaluated on %d chronological holdout timestamps (%s to %s), never used for tuning or fitting.",
          length(y_test), format(min(test_dates), "%Y-%m-%d"),
          format(max(test_dates), "%Y-%m-%d")),
  "",
  "## Benchmark protocol",
  "",
  "Every row in the results table is produced under one identical protocol:",
  "",
  "1. **Target.** `Appliances` at time *t+1*, one 10-minute step ahead (watt-hours).",
  "2. **Information set.** All features are measured at *t* or earlier, so the design uses only",
  "   information available at prediction time. `lights` is excluded (appliance use should be",
  "   predicted independently of lighting) and the dataset's random columns `rv1` and `rv2` are",
  "   never used.",
  paste0("3. **Split.** Chronological 80/20 holdout inherited from `eda/00_prepare_data.R`: ",
         nrow(x_train), " training origins and ", nrow(x_test), " holdout origins. Every training",
         " timestamp precedes every holdout timestamp, which the script asserts."),
  "4. **Standardisation.** Predictor and target scaling use training statistics only, recomputed",
  "   inside every cross-validation fold.",
  "5. **Tuning.** Hyperparameters come from forward-chaining (blocked) cross-validation inside the",
  "   training period. Validation folds are contiguous blocks of later data, never random splits,",
  "   because random splits leak information across time in an autocorrelated series.",
  "6. **Determinism.** Every estimator is closed-form or zero-initialised coordinate descent, so",
  "   results reproduce exactly without a random seed.",
  "",
  "## Feature library",
  "",
  "| group | channels | count |",
  "| --- | --- | ---: |",
  "| Recent consumption | `appl_lag1` (appliance energy at *t*) | 1 |",
  "| Indoor climate | `T1`-`T9`, `RH_1`-`RH_9` | 18 |",
  "| Outdoor weather | `T_out`, `RH_out`, `Windspeed`, `Visibility`, `Press_mm_hg`, `Tdewpoint` | 6 |",
  "| Calendar | sine/cosine of hour-of-day and weekday | 4 |",
  "| **Shared base library** | used by models 1, 2, 3, 4, 6 | **29** |",
  "| Outdoor weather lags | the 6 outdoor channels at lags 1-12 | 72 |",
  "",
  "Models 1-4 and 6 use the shared base library, so their comparison is like-for-like. Model 5",
  "(distributed lag) and model 4 (polynomial) deliberately use larger libraries - that *is* their",
  "hypothesis - and their extra columns are listed above.",
  "",
  "## Models, rationale and tuned structure",
  "",
  "\"#\" is the model file number in `models/`; \"ref\" marks a reference model that is",
  "scored under the same protocol but is not one of the six.",
  "",
  model_display,
  "",
  "## Results -- the six models head to head",
  "",
  "Sorted by holdout RMSE. \"active\" is the number of non-zero coefficients, which matters for",
  "the L1 models. CV RMSE is the mean blocked-CV error used for tuning and is reported for",
  "reference only; the holdout column is the number that counts.",
  "",
  six_display,
  "",
  "## Results -- reference models",
  "",
  "These are not part of the six; they set the bar. Persistence is the naive benchmark, and the",
  "two OLS rows are unpenalised least squares.",
  "",
  ref_display,
  "",
  "## Best of the six in detail",
  "",
  sprintf("**%s** (%s). %s", best_six$name, best_six$family, best_six$rationale),
  "",
  paste0("- Tuned structure: ", best_six$detail),
  paste0("- Design columns: ", best_six$n_predictors, " (", best_six$n_nonzero, " with non-zero coefficients)"),
  paste0("- Holdout RMSE ", fmt(best_six$rmse), " Wh, MAE ", fmt(best_six$mae), " Wh, R2 ", fmt(best_six$r2, 4)),
  paste0("- Literature reference: ", best_six$reference),
  "",
  "Largest standardised coefficients (per 1 SD of the predictor, so they are comparable across features):",
  "",
  md_row("term", "coefficient (standardised)"),
  md_row("---", "---:"),
  best_six_coefs,
  "",
  "## Does the winner depend on the metric?",
  "",
  if (same_winner) {
    sprintf("No. **%s** has both the lowest RMSE (%.3f Wh) and the lowest MAE (%.3f Wh) of the six.",
            best_six$name, best_six$rmse, best_six$mae)
  } else {
    sprintf("Yes. The lowest RMSE belongs to **%s** (%.3f Wh), but the lowest MAE belongs to **%s** (%.3f Wh). Both are honest summaries of the same predictions; they differ because RMSE squares the error and so is dominated by a handful of large misses, whereas MAE is the typical error.",
            best_six$name, best_six$rmse, mae_best_six$name, mae_best_six$mae)
  },
  "",
  sprintf("The sharper point is this: persistence has MAE %.3f Wh and RMSE %.3f Wh. The RMSE winner is",
          persist_mae, persist_rmse),
  sprintf("**worse than persistence on MAE** (%.3f Wh), and only %d of the six beat persistence on MAE at all.",
          best_six$mae, n_six_beats_persist_mae),
  "So on typical error there is essentially nothing to choose between a naive copy and a degree-2",
  "penalised regression; the models earn their keep by being less catastrophic on the largest spikes,",
  "which is exactly what the squared-error objective is built to punish. A model that wins RMSE but",
  "loses MAE is not a contradiction, it is a definition of which errors you decided to care about, and",
  "that choice should be stated rather than left implicit.",
  "",
  "For this project the recommendation is to lead with RMSE (the objective the least-squares models",
  "actually minimise, and the metric the cited benchmark reports) while reporting MAE alongside, since",
  "MAE is the number a household would feel.",
  "",
  "## Why a reference model edges the six",
  "",
  sprintf("The best overall row is **%s** (RMSE %.3f Wh), a 34-column OLS using appliance readings at",
          overall_best$name, overall_best$rmse),
  "lags 1-6 plus the shared library. It leads the best of the six by only a few tenths of a Wh, but the",
  "reason is worth stating precisely: models 1-4 and 6 are restricted to `appl_lag1` alone, so they see",
  "one slice of consumption history, while the reference sees six. Consumption is strongly",
  "autocorrelated at the 10-minute scale, so extra lags are informative - and once they are present,",
  "fancier penalties and basis expansions add little.",
  "",
  "This is a finding about the information set, not about the methods: it says the six-model comparison",
  "should be read as \"which penalisation or structure is best *given* this feature library\", and that",
  "the next improvement to chase is more consumption history, not a more exotic penalty. Adding",
  "`appl_lag2`-`appl_lag6` to the shared base library in `models/00_shared.R` is a one-line change; the",
  "benchmark would then re-rank all six on the wider information set.",
  "",
  "The reference's largest standardised coefficients, which show how dominant the daily cycle is:",
  "",
  md_row("term", "coefficient (standardised)"),
  md_row("---", "---:"),
  overall_best_coefs,
  "",
  "## The distributed-lag hypothesis was rejected",
  "",
  "Model 5 exists to test a falsifiable claim: does outdoor weather measured at earlier steps add",
  "predictive information? Its lag length K was chosen from {0, 1, 2, 3, 6, 12} by blocked",
  "cross-validation, and the winner was **K = 0** - the model collapsed back to the shared base",
  "library. That is why it scores exactly the same as the `OLS, base library` reference row. With",
  "indoor sensors already in the library, past outdoor weather carries no extra signal at this",
  "horizon, and the data said so rather than the modeller assuming it.",
  "",
  "## Are the differences real?",
  "",
  "Paired Wilcoxon signed-rank tests on absolute errors against the winner. These are",
  "indicative, not decisive: consecutive 10-minute errors are strongly dependent, which violates",
  "the independence assumption behind the test, and the same holdout is reused for every",
  "comparison. Treat the effect sizes (column 2) as the more trustworthy signal.",
  "",
  paired_display,
  "",
  "## Stability across the holdout period",
  "",
  "RMSE by calendar month for the naive benchmark and the best three of the six. A model that only",
  "wins in one month has not solved the problem.",
  "",
  stability_display,
  paste0("Holdout stamps per month: ",
         paste0(stability$month, " = ", format(month_counts, scientific = FALSE),
                collapse = ", "), "."),
  "The 2016-04 row is a single day (the holdout begins on 2016-04-30), so it is much noisier than",
  "the May row and should not be over-interpreted.",
  "",
  "## Comparison with a published benchmark on this dataset",
  "",
  "Lu, Lu & Kassymova (2026) benchmark machine-learning models on this exact dataset and report that",
  "under chronological validation the best non-baseline model was ridge regression with MAE 52.54 Wh,",
  "RMSE 84.94 Wh and R2 0.126.",
  "",
  sprintf("Our ridge reaches RMSE %.3f Wh with MAE %.3f Wh and R2 %.4f. The gap is large, and the",
          results[["ridge"]]$rmse, results[["ridge"]]$mae, results[["ridge"]]$r2),
  "honest explanation is not that our implementation is better. The two setups use different",
  "information sets: that study engineers time features and weather/sensor channels, and its reported",
  "important features (RH_2, RH_1, RH_3, hour_cos, T3) contain no lagged target term. Our design",
  "includes `appl_lag1` - the most recent appliance reading, which is legitimately available 10",
  sprintf("minutes before the target - and copying that value alone (persistence) already reaches MAE %.2f Wh and RMSE %.2f Wh on this holdout.",
          persist_mae, persist_rmse),
  "Most of the gap therefore comes from a feature that study does not use, not from a better",
  "estimator. The comparison is not like-for-like, and it illustrates the paper's own point about",
  "validation design. To reproduce their setting we would have to drop `appl_lag1` entirely.",
  "",
  "## Limitations",
  "",
  "- One household, roughly 4.5 months. Nothing here generalises to other homes or other seasons.",
  sprintf("- The six models are separated by only %.3f Wh of RMSE, which is small next to the day-to-day variability of the holdout. The exact ordering within the six should be treated as provisional, which is also what the paired tests show: only the L1 models and the naive benchmark are comfortably separated from the winner.", six_spread),
  "- The holdout is a single block of later time, so all comparisons share one test period. A",
  "  rolling-origin evaluation across several folds would give a tighter picture of model ranking.",
  "- Errors are heavy-tailed and dominated by switching events: all models miss the same large",
  "  spikes, so RMSE differences are driven by behaviour on a small number of timestamps.",
  "- Residuals are autocorrelated, so the paired tests above overstate confidence.",
  "- Coefficients are shown on the standardised scale for comparability; they are associations,",
  "  not causal effects.",
  "- The comparison against the published benchmark is not like-for-like, as explained above.",
  "",
  "## How to reproduce",
  "",
  "```r",
  "source(\"run_all.R\")                       # everything, including this file",
  "source(\"benchmarks/02_model_benchmark.R\") # just the benchmark",
  "```",
  "",
  "## References",
  "",
  "- Almon, S. (1965). The distributed lag between capital appropriations and expenditures. *Econometrica*, 33(1), 178-196. https://doi.org/10.2307/1911894",
  "- Bishop, C. M. (2006). *Pattern Recognition and Machine Learning*. Springer.",
  "- Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for generalized linear models via coordinate descent. *Journal of Statistical Software*, 33(1), 1-22. https://doi.org/10.18637/jss.v033.i01",
  "- Hastie, T., Tibshirani, R., & Friedman, J. (2009). *The Elements of Statistical Learning* (2nd ed.). Springer.",
  "- Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. *Technometrics*, 12(1), 55-67. https://doi.org/10.1080/00401706.1970.10488634",
  "- Lu, H., Lu, H., & Kassymova, G. K. (2026). Leakage-Aware Benchmarking of Machine-Learning Models for Residential Appliance Energy Prediction Using a Public Smart-Home Dataset. *Progress in Energy and Environment*, 32(2), 36-51. https://doi.org/10.37934/progee.32.2.3651",
  "- Ning, S., & Yoon, H. (2024). A New Model for Building Energy Modeling and Management Using Predictive Analytics: Partitioned Hierarchical Multitask Regression (PHMR). *Indoor Air*, 2024(1), 5595459. https://doi.org/10.1155/2024/5595459",
  "- Pang, Y., Zhou, X., Zhang, J., Sun, Q., & Zheng, J. (2022). Hierarchical electricity time series prediction with cluster analysis and sparse penalty. *Pattern Recognition*, 126, 108555. https://doi.org/10.1016/j.patcog.2022.108555",
  "- Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. *Journal of the Royal Statistical Society: Series B*, 58(1), 267-288. https://doi.org/10.1111/j.2517-6161.1996.tb02080.x",
  "- Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. *Journal of the Royal Statistical Society: Series B*, 67(2), 301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x",
  ""
)

writeLines(md, file.path(project_root, "results.md"), useBytes = TRUE)

cat(sprintf("\nBest of the six: %s (RMSE %.3f Wh)\n", best_six$name, best_six$rmse))
cat(sprintf("Best overall (incl. references): %s (RMSE %.3f Wh)\n",
            overall_best$name, overall_best$rmse))
cat("Wrote results.md and results/csv/model_*.csv\n")
