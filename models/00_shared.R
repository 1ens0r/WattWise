# WattWise -- shared library for the six linear models
#
# Provides data loading, leakage-free feature engineering, blocked
# (forward-chaining) cross-validation, the from-scratch linear-algebra solvers
# used by the models, and metric helpers.
#
# Model interface expected by `benchmarks/02_model_benchmark.R`. Each model file
# assigns `wattwise_model <- list(...)` containing:
#   name      : short label
#   family    : method family
#   rationale : one sentence on why the method fits this data
#   reference : literature reference for the method
#   tune      : function(x, y) -> list(hyper = list(...), cv_rmse = numeric)
#   fit       : function(x, y, hyper) -> list(predict = function(xnew) numeric,
#                                              info = list(...))
# `x` is always the candidate matrix from `build_candidates()` and each model
# selects the columns it needs. Tuning and fitting only ever receive training
# rows.
#
# Reproducibility: every estimator here is deterministic (closed-form solves and
# zero-initialised coordinate descent), so no random seed is required.

# ---------------------------------------------------------------------------
# Data and feature engineering
# ---------------------------------------------------------------------------

WATTWISE_INDOOR <- c(paste0("T", 1:9), paste0("RH_", 1:9))
WATTWISE_OUTDOOR <- c("T_out", "RH_out", "Windspeed", "Visibility",
                      "Press_mm_hg", "Tdewpoint")
WATTWISE_ENV_LAG_VARS <- WATTWISE_OUTDOOR

load_series <- function(project_root) {
  path <- file.path(project_root, "data", "processed", "energydata_clean.csv")
  if (!file.exists(path)) stop("Missing ", path, ". Run eda/00_prepare_data.R first.",
                               call. = FALSE)
  d <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  d$date <- as.POSIXct(d$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (anyNA(d$date)) stop("Could not parse timestamps in ", path, call. = FALSE)
  d <- d[order(d$date), , drop = FALSE]
  rownames(d) <- NULL
  d
}

# Column groups. `base` is the shared predictor library; `partition_only`
# columns describe the observation (hour of day, outdoor temperature) and are
# never used as predictors.
wattwise_columns <- function() {
  list(
    indoor = WATTWISE_INDOOR,
    outdoor = WATTWISE_OUTDOOR,
    appliance_lag = paste0("appl_lag", 1:6),
    base = c("appl_lag1", WATTWISE_INDOOR, WATTWISE_OUTDOOR,
             "hour_sin", "hour_cos", "weekday_sin", "weekday_cos"),
    env_lag_vars = WATTWISE_ENV_LAG_VARS,
    partition_only = c("hour", "t_out")
  )
}

wattwise_env_lag_cols <- function(K) {
  if (K <= 0L) return(character(0))
  as.vector(vapply(seq_len(K), function(l) paste0(WATTWISE_ENV_LAG_VARS, "_lag", l),
                   character(length(WATTWISE_ENV_LAG_VARS))))
}

# Build the candidate design matrix.
#
# One row per prediction origin t (target = Appliances at t+1, 10 minutes ahead).
# Every feature is measured at t or earlier, so nothing from the target time is
# used. `lights`, `rv1` and `rv2` are deliberately excluded.
build_candidates <- function(d, env_lags = 12L, appl_lags = 6L) {
  n <- nrow(d)
  max_lag <- max(env_lags, appl_lags)
  origins <- seq.int(max_lag + 1L, n - 1L)
  target_idx <- origins + 1L

  hour <- as.integer(format(d$date[origins], "%H", tz = "UTC"))
  minute <- as.integer(format(d$date[origins], "%M", tz = "UTC"))
  hour_frac <- hour + minute / 60
  weekday <- as.integer(format(d$date[origins], "%u", tz = "UTC"))

  cols <- list()
  for (l in seq_len(appl_lags)) {
    cols[[paste0("appl_lag", l)]] <- d$Appliances[origins - (l - 1L)]
  }
  for (v in WATTWISE_INDOOR) cols[[v]] <- d[[v]][origins]
  for (v in WATTWISE_OUTDOOR) cols[[v]] <- d[[v]][origins]
  cols[["hour_sin"]] <- sin(2 * pi * hour_frac / 24)
  cols[["hour_cos"]] <- cos(2 * pi * hour_frac / 24)
  cols[["weekday_sin"]] <- sin(2 * pi * weekday / 7)
  cols[["weekday_cos"]] <- cos(2 * pi * weekday / 7)
  for (l in seq_len(env_lags)) {
    for (v in WATTWISE_ENV_LAG_VARS) {
      cols[[paste0(v, "_lag", l)]] <- d[[v]][origins - l]
    }
  }
  cols[["hour"]] <- as.numeric(hour)
  cols[["t_out"]] <- d$T_out[origins]

  x <- as.matrix(as.data.frame(cols, check.names = FALSE))
  storage.mode(x) <- "double"

  list(
    x = x,
    y = d$Appliances[target_idx],
    date = d$date[target_idx],
    origin_date = d$date[origins],
    env_lags = env_lags,
    appl_lags = appl_lags,
    max_lag = max_lag
  )
}

# Chronological train/test labelling that reproduces the Phase 1 split.
label_split <- function(candidate, train, test) {
  train_end <- max(train$date)
  test_start <- min(test$date)
  test_end <- max(test$date)

  split <- ifelse(candidate$date >= test_start, "test", "train")

  if (any(candidate$date[split == "test"] < test_start) ||
      any(candidate$date[split == "test"] > test_end)) {
    stop("Test labels fall outside the holdout period.", call. = FALSE)
  }
  if (any(candidate$date[split == "train"] > train_end)) {
    stop("Training labels extend past the training period.", call. = FALSE)
  }
  if (max(candidate$date[split == "train"]) >= min(candidate$date[split == "test"])) {
    stop("Training and test labels overlap.", call. = FALSE)
  }
  split
}

# ---------------------------------------------------------------------------
# Standardisation (training statistics only)
# ---------------------------------------------------------------------------

standardize_fit <- function(x) {
  center <- colMeans(x)
  scale <- apply(x, 2L, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  list(center = center, scale = scale)
}

standardize_apply <- function(x, sc) {
  sweep(sweep(x, 2L, sc$center, "-"), 2L, sc$scale, "/")
}

# Drop zero-variance columns so coordinate descent cannot divide by zero.
drop_constant_columns <- function(x) {
  keep <- apply(x, 2L, stats::sd) > 0
  keep[!is.finite(keep)] <- FALSE
  x[, keep, drop = FALSE]
}

# ---------------------------------------------------------------------------
# Blocked cross-validation (forward chaining, training rows only)
# ---------------------------------------------------------------------------

blocked_blocks <- function(n, n_blocks) {
  unname(split(seq_len(n), cut(seq_len(n), breaks = n_blocks, labels = FALSE)))
}

# Mean validation RMSE over expanding-window folds: fold b trains on blocks
# 1..b-1 and validates on block b. Fold 1 is always training data only, so the
# model never validates on the earliest, least-informed portion.
cv_mean_rmse <- function(x, y, n_blocks, fit_fn, pred_fn) {
  blocks <- blocked_blocks(length(y), n_blocks)
  scores <- numeric(n_blocks - 1L)
  for (b in 2:n_blocks) {
    tr <- unlist(blocks[seq_len(b - 1L)], use.names = FALSE)
    va <- blocks[[b]]
    state <- fit_fn(x[tr, , drop = FALSE], y[tr])
    pred <- pred_fn(state, x[va, , drop = FALSE])
    scores[b - 1L] <- sqrt(mean((y[va] - pred)^2))
  }
  mean(scores)
}

# ---------------------------------------------------------------------------
# Parallel tuning helpers
# ---------------------------------------------------------------------------
# Tuning fits many independent models: one per validation fold, or one per
# (alpha, fold) pair. Those fits share no state and are deterministic, so they
# are the natural place to use the machine's cores. A PSOCK cluster is used
# because the pipeline must also work on Windows, where fork-based `mclapply`
# is unavailable.
#
# Correctness is independent of the number of cores: each task is scored on its
# own and the results are combined by summing in a fixed order, so a run on one
# core and a run on eight produce identical tables. If a cluster cannot be
# created for any reason the work falls back to `lapply` on the current process.

# The environment the library was sourced into. Worker functions are defined at
# the top level of this file, so exporting their names into a worker's global
# environment makes their dependencies resolvable there too.
wattwise_lib_env <- environment()

wattwise_cores <- function() {
  cores <- tryCatch(parallel::detectCores(logical = TRUE),
                    error = function(e) 1L)
  if (length(cores) != 1L || !is.finite(cores) || cores < 2L) return(1L)
  max(1L, min(as.integer(getOption("wattwise.cores", cores - 1L)), 8L))
}

# Map `fn` over `tasks`, optionally in parallel. `fn` and every function it calls
# must be named in `export`; callers pass the constants below.
wattwise_lapply <- function(tasks, fn, export = character(0)) {
  n_tasks <- length(tasks)
  if (n_tasks < 2L) return(lapply(tasks, fn))

  cores <- min(wattwise_cores(), n_tasks)
  if (cores < 2L) return(lapply(tasks, fn))

  out <- tryCatch({
    cluster <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(cluster, varlist = unique(export),
                            envir = wattwise_lib_env)
    parallel::parLapply(cluster, tasks, fn)
  }, error = function(e) NULL)

  if (is.null(out)) lapply(tasks, fn) else out
}

# ---------------------------------------------------------------------------
# Linear-algebra solvers
# ---------------------------------------------------------------------------

soft_threshold <- function(z, gamma) {
  ifelse(z > gamma, z - gamma, ifelse(z < -gamma, z + gamma, 0))
}

# Ridge with an unpenalised intercept, solved on centre-standardised inputs.
ridge_solve <- function(xtx, xty, lambda) {
  as.vector(solve(xtx + diag(lambda, ncol(xtx)), xty))
}

# Ridge penalty grid anchored to the sample size.
#
# The ridge objective is a *sum* of squared errors, so the Gram matrix has
# diagonal about n and the same penalty means something completely different on
# 200 rows and 15000 rows. Expressing the grid as `lambda / n` in
# [lo, hi] makes the range mean the same thing at any sample size: `lambda / n`
# is the fraction of the total variance the penalty removes at full shrinkage,
# so the grid spans "essentially unpenalised" (1e-5) to "almost fully shrunk"
# (10). This is the same normalisation glmnet uses when it stores lambda/n.
ridge_lambda_grid <- function(n_obs, n_lambda = 30L, lo = 1e-5, hi = 10) {
  exp(seq(log(lo * n_obs), log(hi * n_obs), length.out = n_lambda))
}

# One forward-chaining fold of a ridge path. Kept as a top-level function so it
# can be shipped to a worker core by name.
ridge_fold_fit <- function(task) {
  xtr <- drop_constant_columns(task$x_tr)
  sc <- standardize_fit(xtr)
  xt <- standardize_apply(xtr, sc)
  y_mean <- mean(task$y_tr)
  y_sd <- stats::sd(task$y_tr)

  xtx <- crossprod(xt)
  xty <- crossprod(xt, (task$y_tr - y_mean) / y_sd)
  xv <- standardize_apply(task$x_va[, colnames(xtr), drop = FALSE], sc)

  sse <- vapply(task$lambdas, function(l) {
    pred <- as.vector(xv %*% ridge_solve(xtx, xty, l)) * y_sd + y_mean
    sum((task$y_va - pred)^2)
  }, numeric(1))

  list(sse = sse, n = length(task$y_va))
}

WATTWISE_RIDGE_EXPORTS <- c("ridge_fold_fit", "ridge_solve",
                            "drop_constant_columns", "standardize_fit",
                            "standardize_apply")

# Validation RMSE for every lambda on a ridge path, by forward-chaining CV over
# expanding windows. `x` is the design the model actually uses, so the caller
# decides the feature library.
#
# Spinning up a cluster costs about a second on Windows, so it is only worth it
# when the fold fits are genuinely expensive. The dominant cost of a fold is the
# crossproduct of the design, about `nrow * ncol^2` floating-point operations, so
# the work is dispatched across cores only when that exceeds
# `ridge_parallel_flops`. The 29-channel base library is far below the threshold
# and stays serial; the 464-term degree-2 expansion is far above it.
ridge_parallel_flops <- 1e9

ridge_cv_blocks <- function(x, y, lambdas, n_blocks = 6L) {
  blocks <- blocked_blocks(length(y), n_blocks)

  tasks <- lapply(seq_len(n_blocks - 1L), function(i) {
    b <- i + 1L
    tr <- unlist(blocks[seq_len(b - 1L)], use.names = FALSE)
    list(x_tr = x[tr, , drop = FALSE], y_tr = y[tr],
         x_va = x[blocks[[b]], , drop = FALSE], y_va = y[blocks[[b]]],
         lambdas = lambdas)
  })

  folds <- if (nrow(x) * ncol(x)^2 >= ridge_parallel_flops) {
    wattwise_lapply(tasks, ridge_fold_fit, WATTWISE_RIDGE_EXPORTS)
  } else {
    lapply(tasks, ridge_fold_fit)
  }

  sse <- Reduce(`+`, lapply(folds, function(f) f$sse))
  n_obs <- sum(vapply(folds, function(f) f$n, numeric(1)))
  sqrt(sse / n_obs)
}

# Elastic-net path by coordinate descent with warm starts.
#
# minimises (1/(2n)) ||y - x beta||^2 + lambda [ alpha ||beta||_1
#                                                + (1 - alpha)/2 ||beta||_2^2 ]
#
# Reference: Friedman, Hastie & Tibshirani (2010), "Regularization Paths for
# Generalized Linear Models via Coordinate Descent", JSS 33(1).
# Coordinate descent is run to a tolerance of 1e-6 on the largest coefficient
# movement between sweeps, which is far tighter than the noise floor of the
# prediction problem; benchmarks/03_solver_checks.R re-runs it at 1e-10 to
# confirm the fixpoint really satisfies the KKT conditions.
enet_path <- function(x, y, alpha, lambdas, tol = 1e-6, max_iter = 200L) {
  n <- nrow(x)
  p <- ncol(x)
  xx <- colSums(x^2) / n
  beta <- numeric(p)
  residual <- y - as.vector(x %*% beta)
  coefs <- matrix(0, nrow = length(lambdas), ncol = p)

  for (li in seq_along(lambdas)) {
    lam <- lambdas[li]
    for (iter in seq_len(max_iter)) {
      max_delta <- 0
      for (j in seq_len(p)) {
        residual <- residual + x[, j] * beta[j]
        rho <- sum(x[, j] * residual) / n
        new_beta <- soft_threshold(rho, lam * alpha) /
          (xx[j] + lam * (1 - alpha))
        residual <- residual - x[, j] * new_beta
        max_delta <- max(max_delta, abs(new_beta - beta[j]))
        beta[j] <- new_beta
      }
      if (max_delta < tol) break
    }
    coefs[li, ] <- beta
  }
  coefs
}

# Smallest lambda that zeroes every coefficient on a standardised design.
lambda_max <- function(x, y, alpha) {
  max(abs(crossprod(x, y) / nrow(x))) / alpha
}

enet_lambda_grid <- function(x, y, alpha, n_lambda = 30L, ratio = 1e-4) {
  lmax <- lambda_max(x, y, alpha)
  if (!is.finite(lmax) || lmax <= 0) lmax <- 1
  exp(seq(log(lmax), log(lmax * ratio), length.out = n_lambda))
}

# One forward-chaining fold of an elastic-net path, with the fold's own
# standardisation statistics. Top-level so it can be shipped to a worker core.
enet_fold_fit <- function(task) {
  xtr <- drop_constant_columns(task$x_tr)
  sc <- standardize_fit(xtr)
  xt <- standardize_apply(xtr, sc)
  y_mean <- mean(task$y_tr)
  y_sd <- stats::sd(task$y_tr)

  coefs <- enet_path(xt, (task$y_tr - y_mean) / y_sd, task$alpha, task$lambdas)
  xv <- standardize_apply(task$x_va[, colnames(xtr), drop = FALSE], sc)
  pred <- (xv %*% t(coefs)) * y_sd + y_mean

  list(sse = colSums((task$y_va - pred)^2), n = length(task$y_va))
}

WATTWISE_ENET_EXPORTS <- c("enet_fold_fit", "enet_path", "soft_threshold",
                           "drop_constant_columns", "standardize_fit",
                           "standardize_apply")

# Blocked-CV tuning for the elastic-net family. Returns one row per
# (alpha, lambda) with the mean validation RMSE across expanding-window folds.
# The lambda grid is fixed once from the full training set so that every fold is
# scored on the same path, following the usual glmnet CV convention.
#
# Each (alpha, fold) pair is an independent path fit, so they are dispatched
# across the available cores; the validation errors are then summed in a fixed
# order, which makes the returned table independent of the core count.
enet_cv <- function(x, y, alphas, n_blocks = 6L, n_lambda = 25L) {
  blocks <- blocked_blocks(length(y), n_blocks)

  # The lambda grid must be derived from standardised columns, otherwise
  # lambda_max is dominated by whichever channel has the largest raw units.
  x_kept <- drop_constant_columns(x)
  x_std <- standardize_apply(x_kept, standardize_fit(x_kept))
  y_std <- (y - mean(y)) / stats::sd(y)

  keep_cols <- colnames(x_kept)
  grids <- lapply(alphas, function(a) {
    enet_lambda_grid(x_std, y_std, a, n_lambda = n_lambda)
  })
  names(grids) <- as.character(alphas)

  tasks <- list()
  for (a in alphas) {
    lambdas <- grids[[as.character(a)]]
    for (b in 2:length(blocks)) {
      tr <- unlist(blocks[seq_len(b - 1L)], use.names = FALSE)
      va <- blocks[[b]]
      tasks[[length(tasks) + 1L]] <- list(
        alpha = a, lambdas = lambdas,
        x_tr = x[tr, keep_cols, drop = FALSE], y_tr = y[tr],
        x_va = x[va, keep_cols, drop = FALSE], y_va = y[va]
      )
    }
  }

  fits <- wattwise_lapply(tasks, enet_fold_fit, WATTWISE_ENET_EXPORTS)

  # Accumulate in task order, so the sum is bit-for-bit the same as the serial
  # loop it replaces.
  sse <- list()
  n_obs <- list()
  for (i in seq_along(tasks)) {
    key <- as.character(tasks[[i]]$alpha)
    if (is.null(sse[[key]])) {
      sse[[key]] <- numeric(length(tasks[[i]]$lambdas))
      n_obs[[key]] <- 0L
    }
    sse[[key]] <- sse[[key]] + fits[[i]]$sse
    n_obs[[key]] <- n_obs[[key]] + fits[[i]]$n
  }

  out <- lapply(alphas, function(a) {
    key <- as.character(a)
    data.frame(alpha = a, lambda = grids[[key]],
               cv_rmse = sqrt(sse[[key]] / n_obs[[key]]))
  })

  do.call(rbind, out)
}

# OLS via QR, with aliased columns set to zero so predictions stay finite.
ols_fit <- function(x, y) {
  fit <- stats::lm.fit(x, y)
  beta <- fit$coefficients
  beta[is.na(beta)] <- 0
  list(beta = as.vector(beta), aliased = sum(is.na(fit$coefficients)))
}

# Degree-2 polynomial expansion: main effects, squares, pairwise interactions.
expand_poly <- function(x, degree = 2L) {
  p <- ncol(x)
  nm <- colnames(x)
  terms <- list()
  names_out <- character(0)

  add <- function(vec, label) {
    terms[[length(terms) + 1L]] <<- vec
    names_out <<- c(names_out, label)
  }

  for (j in seq_len(p)) add(x[, j], nm[j])
  if (degree >= 2L) {
    for (j in seq_len(p)) add(x[, j]^2, paste0(nm[j], "^2"))
    if (p >= 2L) {
      for (j in seq_len(p - 1L)) {
        for (k in (j + 1L):p) add(x[, j] * x[, k], paste0(nm[j], ":", nm[k]))
      }
    }
  }

  out <- do.call(cbind, terms)
  colnames(out) <- names_out
  out
}

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

rmse_wh <- function(obs, pred) sqrt(mean((obs - pred)^2))
mae_wh <- function(obs, pred) mean(abs(obs - pred))
r2_holdout <- function(obs, pred) {
  1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
}
