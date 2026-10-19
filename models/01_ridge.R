# WattWise -- model 01: Ridge regression (L2-penalised least squares)
#
# Rationale for this dataset: indoor temperature and humidity sensors, and the
# outdoor weather variables, are strongly correlated with each other. Ordinary
# least squares then produces high-variance coefficients, and ridge stabilises
# them by trading a little bias for a large variance reduction. It also gives
# every predictor a non-zero weight, which suits a prediction task where no
# sensor is obviously irrelevant.
#
# Reference: Hoerl, A. E., & Kennard, R. W. (1970). Ridge regression: Biased
#   estimation for nonorthogonal problems. Technometrics, 12(1), 55-67.
#   https://doi.org/10.1080/00401706.1970.10488634
# Dataset-specific evidence: Lu, H., Lu, H., & Kassymova, G. K. (2026).
#   Leakage-Aware Benchmarking of Machine-Learning Models for Residential
#   Appliance Energy Prediction Using a Public Smart-Home Dataset. Progress in
#   Energy and Environment, 32(2), 36-51. https://doi.org/10.37934/progee.32.2.3651
#   -- reports ridge as the best non-baseline model on this dataset under
#   chronological validation (MAE 52.54, RMSE 84.94, R2 0.126).
#
# Implemented with the closed-form ridge solution on centre-standardised
# predictors and a centre-standardised target, so lambda is on a comparable
# scale for every model in the benchmark. Lambda is chosen by blocked
# (forward-chaining) cross-validation on training rows only.

wattwise_model <- list(
  name = "Ridge",
  family = "Linear, L2 penalty",
  rationale = "Stabilises coefficients when weather and indoor sensors are strongly correlated.",
  reference = "Hoerl & Kennard (1970), Technometrics 12(1), 55-67",

  tune = function(x, y) {
    cols <- wattwise_columns()$base
    lambdas <- ridge_lambda_grid(length(y), n_lambda = 30L)
    cv <- ridge_cv_blocks(x[, cols, drop = FALSE], y, lambdas, n_blocks = 6L)
    best <- which.min(cv)
    list(hyper = list(lambda = lambdas[best]), cv_rmse = cv[best])
  },

  fit = function(x, y, hyper) {
    cols <- wattwise_columns()$base
    xtr <- drop_constant_columns(x[, cols, drop = FALSE])
    sc <- standardize_fit(xtr)
    xt <- standardize_apply(xtr, sc)
    y_mean <- mean(y)
    y_sd <- stats::sd(y)

    beta <- ridge_solve(crossprod(xt), crossprod(xt, (y - y_mean) / y_sd),
                        hyper$lambda)

    predict <- function(xnew) {
      xn <- standardize_apply(xnew[, colnames(xtr), drop = FALSE], sc)
      as.vector(xn %*% beta) * y_sd + y_mean
    }

    list(
      predict = predict,
      info = list(
        coefficients = stats::setNames(as.vector(beta), colnames(xtr)),
        n_predictors = ncol(xtr),
        n_nonzero = sum(abs(beta) > 1e-8),
        detail = sprintf("lambda = %.4g (lambda/n = %.3g), 30-value path",
                         hyper$lambda, hyper$lambda / length(y))
      )
    )
  }
)
