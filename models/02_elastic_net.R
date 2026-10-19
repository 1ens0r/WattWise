# WattWise -- model 02: Elastic net regression (L1 + L2 penalty)
#
# Rationale for this dataset: the nine indoor temperature channels (T1-T9) and
# nine indoor humidity channels (RH_1-RH_9) are heavily collinear, and ridge
# alone keeps all of them while lasso tends to pick one arbitrarily. The elastic
# net's grouped-selection behaviour keeps correlated sensor blocks together
# while still shrinking uninformative channels towards zero.
#
# Reference: Zou, H., & Hastie, T. (2005). Regularization and variable selection
#   via the elastic net. Journal of the Royal Statistical Society: Series B,
#   67(2), 301-320. https://doi.org/10.1111/j.1467-9868.2005.00503.x
# Solver reference: Friedman, J., Hastie, T., & Tibshirani, R. (2010).
#   Regularization paths for generalized linear models via coordinate descent.
#   Journal of Statistical Software, 33(1), 1-22. https://doi.org/10.18637/jss.v033.i01
#
# Both the mixing parameter alpha and the penalty lambda are selected by blocked
# cross-validation on training rows only.

wattwise_model <- list(
  name = "Elastic net",
  family = "Linear, L1 + L2 penalty",
  rationale = "Keeps collinear indoor sensor groups together while shrinking uninformative channels.",
  reference = "Zou & Hastie (2005), JRSS-B 67(2), 301-320",

  tune = function(x, y) {
    cols <- wattwise_columns()$base
    xb <- x[, cols, drop = FALSE]
    alphas <- c(0.1, 0.5, 0.9)

    cv <- enet_cv(xb, y, alphas, n_blocks = 6L, n_lambda = 25L)
    best <- which.min(cv$cv_rmse)

    list(
      hyper = list(alpha = cv$alpha[best], lambda = cv$lambda[best]),
      cv_rmse = cv$cv_rmse[best]
    )
  },

  fit = function(x, y, hyper) {
    cols <- wattwise_columns()$base
    xtr <- drop_constant_columns(x[, cols, drop = FALSE])
    sc <- standardize_fit(xtr)
    xt <- standardize_apply(xtr, sc)
    y_mean <- mean(y)
    y_sd <- stats::sd(y)

    coefs <- enet_path(xt, (y - y_mean) / y_sd, hyper$alpha, hyper$lambda)
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
        detail = sprintf("alpha = %.2f, lambda = %.4g, %d/%d active",
                         hyper$alpha, hyper$lambda, sum(abs(beta) > 1e-8),
                         ncol(xtr))
      )
    )
  }
)
