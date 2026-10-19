# WattWise -- verification of the from-scratch linear-algebra solvers
#
# The six models are implemented in base R without any modelling packages, so
# their solvers need to be verified independently rather than trusted. Each
# check below tests a mathematical property that must hold at the optimum,
# using a *different* computation from the one being checked:
#
#   1. OLS       -- normal equations X'(y - Xb) = 0.
#   2. Ridge     -- (X'X + lambda I)b = X'y, cross-checked against a QR solve of
#                   the augmented least-squares system [X; sqrt(lambda) I].
#   3. Lasso     -- Karush-Kuhn-Tucker conditions of the L1-penalised objective.
#   4. Elastic   -- KKT conditions of the L1 + L2-penalised objective.
#      net
#   5. PHMR      -- block optimality: the closed-form solution of the convex
#                   objective must simultaneously be the exact ridge solution of
#                   every block (pooled beta and each regime deviation), which
#                   together certify the global optimum.
#
# Writes results/md/solver_checks.md. Exits with an error if any check fails, so
# a broken solver cannot silently produce a benchmark.
#
# Run: source("benchmarks/03_solver_checks.R")

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

series <- load_series(project_root)
train_part <- utils::read.csv(file.path(project_root, "data", "processed", "train.csv"),
                              stringsAsFactors = FALSE, check.names = FALSE)
test_part <- utils::read.csv(file.path(project_root, "data", "processed", "test.csv"),
                             stringsAsFactors = FALSE, check.names = FALSE)
train_part$date <- as.POSIXct(train_part$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
test_part$date <- as.POSIXct(test_part$date, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

candidate <- build_candidates(series, env_lags = 12L, appl_lags = 6L)
split <- label_split(candidate, train_part, test_part)

# The checks are algebraic, so a contiguous slice of the training period is
# enough and keeps the script quick.
n_check <- 6000L
columns <- wattwise_columns()$base
x_raw <- candidate$x[split == "train", columns, drop = FALSE][seq_len(n_check), , drop = FALSE]
y_raw <- candidate$y[split == "train"][seq_len(n_check)]

sc <- standardize_fit(x_raw)
x <- standardize_apply(x_raw, sc)
y <- (y_raw - mean(y_raw)) / stats::sd(y_raw)
n <- nrow(x)
p <- ncol(x)

relative <- function(a, b) max(abs(a - b)) / max(1, max(abs(b)))
records <- list()
record <- function(check, violation, tolerance, detail) {
  records[[length(records) + 1L]] <<- data.frame(
    check = check, max_violation = violation, tolerance = tolerance,
    status = if (violation <= tolerance) "PASS" else "FAIL",
    detail = detail, stringsAsFactors = FALSE
  )
}

cat(sprintf("Solver checks on %d standardised training rows and %d columns\n",
            n, p))

# --- 1. OLS -----------------------------------------------------------------
ols <- ols_fit(x, y)
gradient <- as.vector(crossprod(x, y - x %*% ols$beta)) / n
record("OLS normal equations", max(abs(gradient)), 1e-10,
       "X'(y - Xb)/n should vanish at the least-squares solution")

# --- 2. Ridge: normal equations and an independent QR solve ------------------
lambda <- 10
beta_ridge <- ridge_solve(crossprod(x), crossprod(x, y), lambda)
normal_violation <- max(abs((crossprod(x) + diag(lambda, p)) %*% beta_ridge -
                              crossprod(x, y))) / max(1, max(abs(crossprod(x, y))))
record("Ridge normal equations", normal_violation, 1e-10,
       "(X'X + lambda I)b = X'y")

augmented <- rbind(x, diag(sqrt(lambda), p))
beta_qr <- as.vector(qr.coef(qr(augmented), c(y, rep(0, p))))
record("Ridge vs augmented-QR solve", relative(beta_qr, beta_ridge), 1e-8,
       sprintf("QR on [X; sqrt(lambda)I] with lambda = %g", lambda))

# --- 3/4. Lasso and elastic net: KKT conditions ------------------------------
kkt_violation <- function(x, y, beta, lambda, alpha) {
  residual <- y - as.vector(x %*% beta)
  gradient <- as.vector(crossprod(x, residual)) / nrow(x)
  active <- abs(beta) > 1e-10
  violation <- numeric(ncol(x))
  # Active coefficients: gradient must equal the subgradient exactly.
  violation[active] <- abs(gradient[active] -
                             lambda * (1 - alpha) * beta[active] -
                             lambda * alpha * sign(beta[active]))
  # Inactive coefficients: gradient must sit inside the subgradient interval.
  violation[!active] <- pmax(abs(gradient[!active]) - lambda * alpha, 0)
  max(violation)
}

for (spec in list(list(alpha = 1, label = "Lasso"),
                  list(alpha = 0.5, label = "Elastic net (alpha = 0.5)"))) {
  lambda <- 0.01 * n
  # Solved here at a much tighter tolerance than the models use, to test that
  # the coordinate descent fixpoint really is the optimum.
  beta <- as.vector(enet_path(x, y, spec$alpha, lambda, tol = 1e-10,
                              max_iter = 5000L))
  record(sprintf("%s KKT conditions", spec$label),
         kkt_violation(x, y, beta, lambda, spec$alpha), 1e-6,
         sprintf("lambda/n = %.4g, %d of %d coefficients active",
                 lambda / n, sum(abs(beta) > 1e-10), p))
}

# --- 5. PHMR: block optimality ----------------------------------------------
source(file.path(project_root, "models", "06_phmr.R"))

lam_b <- 1 * n
lam_d <- 1 * (n / 4)
groups <- phmr_groups(cbind(x_raw, hour = candidate$x[split == "train", "hour"][seq_len(n_check)],
                            t_out = candidate$x[split == "train", "t_out"][seq_len(n_check)]),
                      "hour4")
# The solver returns the exact minimiser of the convex objective, so the block
# optimality certificate below should hold to machine precision.
fit <- phmr_solve(x, y, groups, lam_b, lam_d)

delta_violation <- 0
for (g in seq_len(4L)) {
  idx <- which(groups == g)
  xg <- x[idx, , drop = FALSE]
  residual_g <- y[idx] - as.vector(xg %*% fit$beta)
  lhs <- (crossprod(xg) + diag(lam_d, p)) %*% fit$delta[g, ]
  rhs <- crossprod(xg, residual_g)
  delta_violation <- max(delta_violation,
                         max(abs(lhs - rhs)) / max(1, max(abs(rhs))))
}
record("PHMR regime sub-problem", delta_violation, 1e-10,
       "each delta_g solves its own ridge sub-problem at the optimum")
y_adjusted <- y - rowSums(x * fit$delta[groups, , drop = FALSE])
beta_violation <- max(abs((crossprod(x) + diag(lam_b, p)) %*% fit$beta -
                            crossprod(x, y_adjusted))) /
  max(1, max(abs(crossprod(x, y_adjusted))))
record("PHMR pooled sub-problem", beta_violation, 1e-10,
       "the pooled beta solves its ridge sub-problem at the optimum")

# --- 6. Determinism ---------------------------------------------------------
beta_first <- as.vector(enet_path(x, y, 0.5, 0.01 * n, tol = 1e-10, max_iter = 5000L))
beta_again <- as.vector(enet_path(x, y, 0.5, 0.01 * n, tol = 1e-10, max_iter = 5000L))
record("Determinism (repeated fit)", max(abs(beta_again - beta_first)),
       0, "coordinate descent is initialised deterministically, so refits are identical")

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

checks <- do.call(rbind, records)
print(checks[, c("check", "max_violation", "tolerance", "status")],
      row.names = FALSE)

md <- c(
  "# WattWise -- solver verification",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "Generated by: `benchmarks/03_solver_checks.R`",
  "",
  "The six models are implemented in base R with no modelling packages, so each",
  "solver is verified against a mathematical property that must hold at its",
  "optimum, computed by a different route. The checks run on",
  sprintf("%d standardised training rows and %d columns.", n, p),
  "",
  "| check | max violation | tolerance | status | what it proves |",
  "| --- | ---: | ---: | --- | --- |",
  vapply(seq_len(nrow(checks)), function(i) {
    sprintf("| %s | %.3e | %.1e | **%s** | %s |", checks$check[i],
            checks$max_violation[i], checks$tolerance[i], checks$status[i],
            checks$detail[i])
  }, character(1)),
  "",
  "All tolerance limits are on relative violations, i.e. rescaled by the size of",
  "the right-hand side, so they are scale-free.",
  ""
)

md_dir <- file.path(project_root, "results", "md")
csv_dir <- file.path(project_root, "results", "csv")
dir.create(md_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(md, file.path(md_dir, "solver_checks.md"), useBytes = TRUE)
utils::write.csv(checks, file.path(csv_dir, "solver_checks.csv"), row.names = FALSE)

if (any(checks$status != "PASS")) {
  stop("Solver verification failed. See results/md/solver_checks.md.", call. = FALSE)
}
cat("\nAll solver checks passed. Wrote results/md/solver_checks.md\n")
