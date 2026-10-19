# WattWise -- model 06: Partitioned hierarchical multitask regression (PHMR)
#
# Rationale for this dataset: one global linear relationship is a strong
# assumption for a house whose behaviour changes through the day. PHMR
# partitions the observations into regimes and estimates a *related* regression
# for each regime, sharing strength between them instead of fitting each subset
# independently. Formally the coefficients are
#
#     beta_g = beta + delta_g ,   g = 1..G
#
# where beta is the pooled coefficient vector shared by every regime and delta_g
# is a regime-specific deviation shrunk towards zero. The objective is
#
#     min sum_g ||y_g - X_g (beta + delta_g)||^2
#         + lambda_b ||beta||^2 + lambda_d sum_g ||delta_g||^2
#
# so large lambda_b forces one global model and large lambda_d forces full
# pooling, while intermediate values let regimes deviate only where the data
# justify it.
#
# The objective is a convex quadratic in the (G + 1) * p = 145 unknown
# coefficients, so instead of iterating it is assembled and solved once in
# closed form: the joint normal equations
#
#     (Z'Z + diag(lambda_b I_p, lambda_d I_p, ...)) eta = Z'y
#
# with Z = [X, X_1, ..., X_G] and eta = (beta, delta_1, ..., delta_G). The block
# structure of Z'Z (the X_g have disjoint row supports, so X_g'X_h = 0 for
# g != h) makes it cheap to assemble. This returns the exact global optimum;
# benchmarks/03_solver_checks.R certifies it against the block optimality
# conditions.
#
# Reference: Ning, S., & Yoon, H. (2024). A New Model for Building Energy
#   Modeling and Management Using Predictive Analytics: Partitioned Hierarchical
#   Multitask Regression (PHMR). Indoor Air, 2024(1), 5595459.
#   https://doi.org/10.1155/2024/5595459
#   (The paper combines recursive partition regression with hierarchical
#   multitask learning; it reports 32.9-41.8% accuracy gains over
#   single-predictor regression in simulation and a 37% gain on a real HVAC case.)
#
# Implementation scope, stated plainly: this is a two-level hierarchical
# (partial-pooling) linear model with a partitioned coefficient structure. It
# follows the shared-plus-deviation idea of PHMR, but it is NOT a reimplementation
# of the paper's full recursive-partition + multitask algorithm. Both the
# partition definition and the two penalties are chosen by blocked
# cross-validation on training rows only.

PHMR_PARTITIONS <- c("hour4", "tout4")

phmr_groups <- function(x, definition, breaks = NULL) {
  if (identical(definition, "hour4")) {
    # 00-05, 06-11, 12-17, 18-23
    return(findInterval(x[, "hour"], c(6, 12, 18)) + 1L)
  }
  if (identical(definition, "tout4")) {
    # outdoor-temperature quartiles estimated on training data
    return(findInterval(x[, "t_out"], breaks) + 1L)
  }
  stop("Unknown PHMR partition definition: ", definition, call. = FALSE)
}

phmr_train <- function(x, y, definition, lam_b, lam_d) {
  cols <- wattwise_columns()$base
  xtr <- drop_constant_columns(x[, cols, drop = FALSE])
  sc <- standardize_fit(xtr)
  xt <- standardize_apply(xtr, sc)
  y_mean <- mean(y)
  y_sd <- stats::sd(y)

  breaks <- if (identical(definition, "tout4")) {
    as.vector(stats::quantile(x[, "t_out"], c(0.25, 0.5, 0.75), names = FALSE))
  } else {
    NULL
  }

  groups <- phmr_groups(x, definition, breaks)
  fit <- phmr_solve(xt, (y - y_mean) / y_sd, groups, lam_b, lam_d)

  list(
    cols = colnames(xtr), sc = sc, y_mean = y_mean, y_sd = y_sd,
    definition = definition, breaks = breaks,
    beta = fit$beta, delta = fit$delta, objective = fit$objective
  )
}

phmr_predict <- function(state, xnew) {
  xn <- standardize_apply(xnew[, state$cols, drop = FALSE], state$sc)
  g <- phmr_groups(xnew, state$definition, state$breaks)
  coef <- state$beta + state$delta[g, , drop = FALSE]
  as.vector(rowSums(xn * coef)) * state$y_sd + state$y_mean
}

# Exact closed-form solve of the joint ridge system.
#
# Returns the global minimiser of
#   ||y - X beta - sum_g X_g delta_g||^2 + lambda_b ||beta||^2
#     + lambda_d sum_g ||delta_g||^2
# by assembling the (G+1)p normal equations. `groups` assigns each observation
# to exactly one regime.
phmr_solve <- function(x, y, groups, lam_b, lam_d) {
  n_groups <- 4L
  p <- ncol(x)
  size <- (n_groups + 1L) * p

  design <- matrix(0, size, size)
  target <- numeric(size)
  beta_idx <- seq_len(p)

  target[beta_idx] <- crossprod(x, y)
  for (g in seq_len(n_groups)) {
    idx <- which(groups == g)
    block <- p * g + seq_len(p)
    xg <- x[idx, , drop = FALSE]
    gram <- crossprod(xg)

    # X'X = sum_g X_g'X_g because the group supports partition the rows.
    design[beta_idx, beta_idx] <- design[beta_idx, beta_idx] + gram
    design[beta_idx, block] <- gram
    design[block, beta_idx] <- gram
    design[block, block] <- gram
    target[block] <- crossprod(xg, y[idx])
  }

  diag(design)[beta_idx] <- diag(design)[beta_idx] + lam_b
  for (g in seq_len(n_groups)) {
    block <- p * g + seq_len(p)
    diag(design)[block] <- diag(design)[block] + lam_d
  }

  eta <- solve(design, target)
  beta <- as.vector(eta[beta_idx])
  delta <- matrix(0, n_groups, p)
  for (g in seq_len(n_groups)) {
    delta[g, ] <- eta[p * g + seq_len(p)]
  }

  residual <- y - as.vector(x %*% beta) -
    rowSums(x * delta[groups, , drop = FALSE])
  objective <- sum(residual^2) + lam_b * sum(beta^2) + lam_d * sum(delta^2)

  list(beta = beta, delta = delta, objective = objective)
}

wattwise_model <- list(
  name = "PHMR",
  family = "Linear, partitioned hierarchical (partial pooling)",
  rationale = "Lets the linear relationship differ by regime while pooling strength across regimes.",
  reference = "Ning & Yoon (2024), Indoor Air 2024(1), 5595459",

  tune = function(x, y) {
    # Penalties are anchored to the sample size because the objective is a sum
    # of squared errors: the Gram matrix has diagonal about n, and a group's
    # Gram matrix about n/G. Expressed this way, 1e-4 is almost no shrinkage and
    # 10 is almost complete shrinkage, at any sample size.
    n <- length(y)
    lam_b_grid <- n * c(1e-4, 1e-2, 1e-1, 1, 10)
    lam_d_grid <- (n / 4) * c(1e-4, 1e-2, 1e-1, 1, 10)
    best <- list(cv_rmse = Inf, definition = NA_character_, lam_b = NA_real_,
                 lam_d = NA_real_)

    for (definition in PHMR_PARTITIONS) {
      for (lam_b in lam_b_grid) {
        for (lam_d in lam_d_grid) {
          cv <- cv_mean_rmse(
            x, y, n_blocks = 5L,
            fit_fn = function(xt, yt) phmr_train(xt, yt, definition, lam_b, lam_d),
            pred_fn = function(state, xv) phmr_predict(state, xv)
          )
          if (cv < best$cv_rmse) {
            best <- list(cv_rmse = cv, definition = definition,
                         lam_b = lam_b, lam_d = lam_d)
          }
        }
      }
    }

    list(
      hyper = list(definition = best$definition, lam_b = best$lam_b,
                   lam_d = best$lam_d),
      cv_rmse = best$cv_rmse
    )
  },

  fit = function(x, y, hyper) {
    state <- phmr_train(x, y, hyper$definition, hyper$lam_b, hyper$lam_d)

    predict <- function(xnew) phmr_predict(state, xnew)

    regime_label <- if (identical(hyper$definition, "hour4")) {
      "hour-of-day blocks 00-05/06-11/12-17/18-23"
    } else {
      sprintf("outdoor-temperature quartiles [%.1f, %.1f, %.1f]",
              state$breaks[1L], state$breaks[2L], state$breaks[3L])
    }

    list(
      predict = predict,
      info = list(
        coefficients = stats::setNames(state$beta, state$cols),
        n_predictors = length(state$cols),
        n_nonzero = sum(abs(state$beta) > 1e-8),
        detail = sprintf("G = 4 regimes by %s, lambda_b/n = %.3g, lambda_d/(n/G) = %.3g, exact solve",
                         regime_label, hyper$lam_b / length(y),
                         hyper$lam_d / (length(y) / 4)),
        extras = list(
          regime_coefficients = state$beta + state$delta,
          regime_labels = regime_label
        )
      )
    )
  }
)
