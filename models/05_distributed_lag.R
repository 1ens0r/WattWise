# WattWise -- model 05: Distributed-lag linear regression
#
# Rationale for this dataset: thermal inertia means the effect of outdoor
# conditions on appliance use is not instantaneous. A distributed-lag model adds
# outdoor weather measured at t-1, t-2, ... as separate regressors, so the
# regression estimates the whole lag profile rather than assuming today's
# weather acts immediately. It tests a falsifiable claim: does earlier
# environmental history add predictive information beyond the current reading?
#
# References:
#   Almon, S. (1965). The distributed lag between capital appropriations and
#     expenditures. Econometrica, 33(1), 178-196.
#     https://doi.org/10.2307/1911894
#   Related energy time-series work: Pang, Y., Zhou, X., Zhang, J., Sun, Q., &
#     Zheng, J. (2022). Hierarchical electricity time series prediction with
#     cluster analysis and sparse penalty. Pattern Recognition, 126, 108555.
#     https://doi.org/10.1016/j.patcog.2022.108555
#
# The lag length K is selected by blocked cross-validation on training rows only.
# K = 0 is in the grid, so the data can reject the lag structure and collapse the
# model back to the shared base library.

wattwise_model <- list(
  name = "Distributed lag (OLS)",
  family = "Linear, unrestricted lag profile",
  rationale = "Tests whether earlier outdoor weather adds information beyond the current reading.",
  reference = "Almon (1965), Econometrica 33(1), 178-196",

  tune = function(x, y) {
    base <- wattwise_columns()$base
    Ks <- c(0L, 1L, 2L, 3L, 6L, 12L)
    scores <- numeric(length(Ks))

    for (i in seq_along(Ks)) {
      cols <- c(base, wattwise_env_lag_cols(Ks[i]))
      xk <- x[, cols, drop = FALSE]
      scores[i] <- cv_mean_rmse(
        xk, y, n_blocks = 6L,
        fit_fn = function(xt, yt) ols_fit(xt, yt),
        pred_fn = function(state, xv) as.vector(xv %*% state$beta)
      )
    }

    best <- which.min(scores)
    list(
      hyper = list(K = Ks[best],
                   columns = c(base, wattwise_env_lag_cols(Ks[best]))),
      cv_rmse = scores[best]
    )
  },

  fit = function(x, y, hyper) {
    xtr <- x[, hyper$columns, drop = FALSE]
    state <- ols_fit(xtr, y)

    predict <- function(xnew) {
      as.vector(xnew[, hyper$columns, drop = FALSE] %*% state$beta)
    }

    list(
      predict = predict,
      info = list(
        coefficients = stats::setNames(state$beta, hyper$columns),
        n_predictors = ncol(xtr),
        n_nonzero = sum(abs(state$beta) > 1e-8),
        detail = sprintf("%d lag step(s) of outdoor weather, %d predictors, %d aliased",
                         hyper$K, ncol(xtr), state$aliased)
      )
    )
  }
)
