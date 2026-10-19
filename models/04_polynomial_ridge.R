# WattWise -- model 04: Polynomial regression with ridge regularisation
#
# Rationale for this dataset: appliance consumption responds to weather and
# indoor conditions in a non-linear way (heating and cooling thresholds, and
# humidity effects that matter only at high temperature). A degree-2 expansion
# adds squared terms and pairwise interactions, which lets the model bend while
# remaining *linear in its coefficients*, so it is still estimated by penalised
# least squares and is still a linear model in the statistical sense. The ridge
# penalty is essential here: the expansion takes the 29-channel library to 464
# terms whose columns are strongly collinear.
#
# References:
#   Hastie, T., Tibshirani, R., & Friedman, J. (2009). The Elements of
#     Statistical Learning (2nd ed.), Ch. 3.1 and 3.4. Springer.
#   Bishop, C. M. (2006). Pattern Recognition and Machine Learning, Sec. 3.1.
#     Springer.
#
# Note: the reference supplied with the project brief (ISBN 978-1-315-37027-4)
# could not be resolved to a verifiable catalogue record, so the two standard
# texts above are cited instead.
#
# Lambda is selected by blocked cross-validation on training rows only.

wattwise_model <- list(
  name = "Polynomial ridge",
  family = "Linear in coefficients, degree-2 expansion",
  rationale = "Captures threshold and interaction effects in weather and indoor sensors under an L2 penalty.",
  reference = "Hastie, Tibshirani & Friedman (2009), ESL 2nd ed., Ch. 3.1/3.4",

  tune = function(x, y) {
    cols <- wattwise_columns()$base
    xp <- expand_poly(x[, cols, drop = FALSE], degree = 2L)
    lambdas <- ridge_lambda_grid(length(y), n_lambda = 30L)
    cv <- ridge_cv_blocks(xp, y, lambdas, n_blocks = 5L)
    best <- which.min(cv)
    list(hyper = list(lambda = lambdas[best]), cv_rmse = cv[best])
  },

  fit = function(x, y, hyper) {
    cols <- wattwise_columns()$base
    xp <- expand_poly(x[, cols, drop = FALSE], degree = 2L)
    xtr <- drop_constant_columns(xp)
    sc <- standardize_fit(xtr)
    xt <- standardize_apply(xtr, sc)
    y_mean <- mean(y)
    y_sd <- stats::sd(y)

    beta <- ridge_solve(crossprod(xt), crossprod(xt, (y - y_mean) / y_sd),
                        hyper$lambda)

    predict <- function(xnew) {
      xpn <- expand_poly(xnew[, cols, drop = FALSE], degree = 2L)
      xn <- standardize_apply(xpn[, colnames(xtr), drop = FALSE], sc)
      as.vector(xn %*% beta) * y_sd + y_mean
    }

    list(
      predict = predict,
      info = list(
        coefficients = stats::setNames(as.vector(beta), colnames(xtr)),
        n_predictors = ncol(xtr),
        n_nonzero = sum(abs(beta) > 1e-8),
        detail = sprintf("degree 2, %d terms (%d main effects), lambda/n = %.3g",
                         ncol(xtr), length(cols), hyper$lambda / length(y))
      )
    )
  }
)
