# WattWise -- model 03: Lasso regression (L1 penalty)
#
# Rationale for this dataset: with 28 sensor and calendar channels, a sparse
# model is useful both for prediction and for interpretation, because it states
# explicitly which measurements matter. The lasso also acts as an embedded
# feature selector, giving a defensible shortlist of predictors rather than
# relying on univariate screening (which would ignore the correlation
# structure).
#
# Reference: Tibshirani, R. (1996). Regression shrinkage and selection via the
#   lasso. Journal of the Royal Statistical Society: Series B, 58(1), 267-288.
#   https://doi.org/10.1111/j.2517-6161.1996.tb02080.x
# Solver reference: Friedman, J., Hastie, T., & Tibshirani, R. (2010).
#   Regularization paths for generalized linear models via coordinate descent.
#   Journal of Statistical Software, 33(1), 1-22. https://doi.org/10.18637/jss.v033.i01
#
# Lambda is selected by blocked cross-validation on training rows only. The
# number of non-zero coefficients is reported in the benchmark table, because it
# is part of the model's claim, not just its accuracy.

wattwise_model <- list(
  name = "Lasso",
  family = "Linear, L1 penalty",
  rationale = "Selects a compact, interpretable subset of the 28 sensor and calendar channels.",
  reference = "Tibshirani (1996), JRSS-B 58(1), 267-288",

  tune = function(x, y) {
    cols <- wattwise_columns()$base
    xb <- x[, cols, drop = FALSE]

    cv <- enet_cv(xb, y, alphas = 1, n_blocks = 6L, n_lambda = 25L)
    best <- which.min(cv$cv_rmse)

    list(hyper = list(lambda = cv$lambda[best]), cv_rmse = cv$cv_rmse[best])
  },

  fit = function(x, y, hyper) {
    cols <- wattwise_columns()$base
    xtr <- drop_constant_columns(x[, cols, drop = FALSE])
    sc <- standardize_fit(xtr)
    xt <- standardize_apply(xtr, sc)
    y_mean <- mean(y)
    y_sd <- stats::sd(y)

    coefs <- enet_path(xt, (y - y_mean) / y_sd, alpha = 1, hyper$lambda)
    beta <- as.vector(coefs[1L, ])

    predict <- function(xnew) {
      xn <- standardize_apply(xnew[, colnames(xtr), drop = FALSE], sc)
      as.vector(xn %*% beta) * y_sd + y_mean
    }

    list(
      predict = predict,
      info = list(
        coefficients = stats::setNames(beta, colnames(xtr)),
        n_predictors = ncol(xtr),
        n_nonzero = sum(abs(beta) > 1e-8),
        detail = sprintf("lambda = %.4g, %d/%d active", hyper$lambda,
                         sum(abs(beta) > 1e-8), ncol(xtr))
      )
    )
  }
)
