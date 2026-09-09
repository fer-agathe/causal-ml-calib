# ______________________________________________________________________________
# Causal Estimators (ATE & GATE)
# ______________________________________________________________________________

library(tidyverse)
library(sandwich)
library(ebal)



#' Simple Regression-Based ATE Estimator
#' @keywords internal
compute_ate_reg <- function(g0, g1) {
  mean(g1 - g0)
}

# Unnormalized IPW and AIPW----

#' Inverse Probability Weighting (IPW) ATE Estimator
#' @keywords internal
compute_ate_ipw <- function(m, d, y) {
  mean( (d * y / m) - ((1 - d) * y / (1 - m)) )
}

#' Doubly Robust AIPW ATE Estimator
#' @keywords internal
compute_ate_aipw <- function(g0, g1, m, d, y) {
  mean( (g1 - g0) + (d * (y - g1) / m) - ((1 - d) * (y - g0) / (1 - m)) )
}

# Hajeck Normalized IPW and AIPW----

#' Normalized / Stabilized Inverse Probability Weighting (Hajek ATE)
#'
#' @param m Numeric vector of predicted propensity scores.
#' @param d Numeric vector of binary treatment indicators.
#' @param y Numeric vector of observed outcomes.
#'
#' @return Numeric scalar of normalized IPW ATE.
#' @export
compute_ate_ipw_normalized <- function(m, d, y) {
  w1 <- d / m
  w0 <- (1 - d) / (1 - m)

  e1 <- sum(w1 * y) / sum(w1)
  e0 <- sum(w0 * y) / sum(w0)

  e1 - e0
}

#' Normalized Doubly Robust / AIPW ATE Estimator
#'
#' @export
compute_ate_aipw_normalized <- function(g0, g1, m, d, y) {
  w1 <- d / m
  w0 <- (1 - d) / (1 - m)

  # Standard baseline outcome regressions
  dr_base <- mean(g1 - g0)

  # Normalized residual corrections
  adj1 <- sum(w1 * (y - g1)) / sum(w1)
  adj0 <- sum(w0 * (y - g0)) / sum(w0)

  dr_base + adj1 - adj0
}

# Equilibrium / Balancing IPW and AIPW----

compute_ate_eqb_reg <- function(y, d, geqb_pred) {
  n0 <- sum(d == 0)
  n1 <- sum(d == 1)

  robust_signal <- (1/n1) * (y - geqb_pred) - (1/n0) * (y - geqb_pred) * (1 - d)
  mean(robust_signal)
}

compute_ate_eqb_ipw <- function(y, d, m_pred) {
  n0 <- sum(d == 0)
  n1 <- sum(d == 1)

  robust_signal <- (1/n1 + 1/n0) * y * (d - m_pred)
  mean(robust_signal)
}

compute_ate_eqb_aipw <- function(y, d, geqb_pred, m_pred) {
  n0 <- sum(d == 0)
  n1 <- sum(d == 1)

  robust_signal <- (1/n1 + 1/n0) * (y - geqb_pred) * (d - m_pred)
  mean(robust_signal)
}

# Fold Aggregator----

#' Aggregate ATE estimates across folds (cross-fitting)
#'
#' @param pred_folds List of fold objects (see [nuisance_estim()]).
#' @param n_folds Number of folds used in the cross-fitting procedure.
#'
#' @description
#' Computes fold-specific ATE estimates using outcome regression, unnormalized
#' and Hajek-normalized IPW/AIPW, and equilibrium estimators, then aggregates
#' them by averaging across folds.
#'
#' @md
compute_ate_folds <- function(pred_folds, n_folds) {
  ate_folds <- vector(mode = "list", length = n_folds)

  for (k in seq_len(n_folds)) {
    g0_pred_k   <- pred_folds[[k]]$g0_pred_k
    g1_pred_k   <- pred_folds[[k]]$g1_pred_k
    geqb_pred_k <- pred_folds[[k]]$geqb_pred_k
    m_pred_k    <- pred_folds[[k]]$m_pred_k
    y_k         <- pred_folds[[k]]$y_k
    d_k         <- pred_folds[[k]]$d_k

    # 1. Baseline Regression Estimator
    ate_reg_k <- compute_ate_reg(g0 = g0_pred_k, g1 = g1_pred_k)

    # 2. Unnormalized Estimators
    ate_ipw_k  <- compute_ate_ipw(m = m_pred_k, d = d_k, y = y_k)
    ate_aipw_k <- compute_ate_aipw(
      g0 = g0_pred_k, g1 = g1_pred_k, m = m_pred_k, d = d_k, y = y_k
    )

    # 3. Hajek Normalized Estimators
    ate_ipw_norm_k <- compute_ate_ipw_normalized(
      m = m_pred_k, d = d_k, y = y_k
    )
    ate_aipw_norm_k <- compute_ate_aipw_normalized(
      g0 = g0_pred_k, g1 = g1_pred_k, m = m_pred_k, d = d_k, y = y_k
    )

    # 4. Equilibrium Estimators (Flachaire et al., 2025)
    ate_eqb_reg_k  <- compute_ate_eqb_reg(
      y = y_k, d = d_k, geqb_pred = geqb_pred_k
    )
    ate_eqb_ipw_k  <- compute_ate_eqb_ipw(
      y = y_k, d = d_k, m_pred = m_pred_k
    )
    ate_eqb_aipw_k <- compute_ate_eqb_aipw(
      y = y_k, d = d_k, geqb_pred = geqb_pred_k, m_pred = m_pred_k
    )

    # Store fold results
    ate_folds[[k]] <- list(
      ate_reg       = ate_reg_k,
      ate_ipw       = ate_ipw_k,
      ate_ipw_norm  = ate_ipw_norm_k,
      ate_aipw      = ate_aipw_k,
      ate_aipw_norm = ate_aipw_norm_k,
      ate_eqb_reg   = ate_eqb_reg_k,
      ate_eqb_ipw   = ate_eqb_ipw_k,
      ate_eqb_aipw  = ate_eqb_aipw_k
    )
  }

  # Average across folds
  ate <- names(ate_folds[[1]]) |>
    map(\(ate_type) mean(map_dbl(ate_folds, ate_type)))
  names(ate) <- names(ate_folds[[1]])

  as_tibble(ate)
}

# GATE----

#' Max-|Z| Quantile for Simultaneous Confidence Intervals
#'
#' Computes the \eqn{(1 - \alpha)} quantile of the distribution of
#' \eqn{\max_{j=1,\dots,p} |Z_j|}, where \eqn{Z_j \sim \mathcal{N}(0,1)}
#' are independent.
#' @param p Dimension of the parameter vector (number of simultaneous statistics).
#' @param n_sim Number of Monte Carlo simulations (default to 10,000).
#' @param alpha Significance level, in (0,1).
#'
#' @returns A single numeric value corresponding to the
#'   \eqn{(1 - \alpha)} quantile of the maximum absolute standard normal
#'   distribution.
qtmax <- function(p,
                  n_sim = 10000,
                  alpha) {
  z <- matrix(rnorm(p * n_sim), nrow = p, ncol = n_sim)
  tmaxs <- apply(abs(z), 2, max)

  quantile(tmaxs, probs = 1 - alpha, names = FALSE)
}

#' Compute Group Average Treatment Effects (GATE) via Quantile Binning + OLS
#'
#' Computes group-specific average effects by partitioning a scalar covariate
#' `X` into quantile-defined bins, regressing `Y` on bin indicators
#' (no intercept), and forming pointwise and simultaneous confidence intervals
#' using robust (HC/White) standard errors.
#'
#' @param X Scalar running/heterogeneity variable used to define quantile bins
#' @param Y Outcome variable to be regressed on the bin indicators.
#' @param max_grid Integer >= 2. Number of quantile segments used to partition `X`.
#' @param alpha ignificance level for confidence intervals, in (0,1). Default to
#' `0.05` corresponds to 95% intervals.
#' @param n_sim Number of Monte Carlo draws used by `qtmax()` to approximate
#' the critical value for the simultaneous confidence band.
#'
#' @returns A list with:
#'   * `beta_hat`: Vector of estimated bin effects (length `max_grid`).
#'   * `ghat_lower_point`, `ghat_upper_point`: Pointwise CI bounds.
#'   * `ghat_lower`, `ghat_upper`: Simultaneous CI bounds.
#'   * `crit_val`: Simultaneous critical value used.
#'   * `grid`: Quantile cutpoints used to form bins (length `max_grid + 1`).
group_average_treatment_effect <- function(X,
                                           Y,
                                           max_grid = 5,
                                           alpha = 0.05,
                                           n_sim = 10000) {
  grid <- quantile(X, probs = c((0:max_grid) / max_grid))
  Xraw <- matrix(NA, nrow = length(Y), ncol = length(grid) - 1)

  # Matrix with quantile-defined bin in which each obs belongs to
  for (k in 2:((length(grid)))) {
    Xraw[, k - 1] <- sapply(X, function(x) ifelse(x >= grid[k - 1] & x < grid[k], 1, 0))
  }
  k <- length(grid)
  Xraw[, k - 1] <- sapply(X, function(x) ifelse(x >= grid[k - 1] & x <= grid[k], 1, 0))

  ols_fit <- lm(Y ~ Xraw - 1)
  coefs <- coef(ols_fit)
  vars <- names(coefs)
  hcv_coefs <- sandwich::vcovHC(ols_fit, type = "HC")
  coefs_se <- sqrt(diag(hcv_coefs)) # White std errors
  ## this is an identity matrix
  c_coefs <- (diag(1 / coefs_se)) %*% hcv_coefs %*% (diag(1 / coefs_se))

  # regression estimate: GATE
  tes <- coefs
  tes_se <- coefs_se
  tes_cor <- c_coefs

  # Simultaneous Confidence Interval
  crit_val <- qtmax(p = ncol(tes_cor), n_sim = n_sim, alpha = 0.05)
  tes_ucb <- tes + crit_val * tes_se
  tes_lcb <- tes - crit_val * tes_se

  # Pointwise Confidence Interval
  tes_uci <- tes + qnorm(1 - alpha / 2) * tes_se
  tes_lci <- tes + qnorm(alpha / 2) * tes_se

  list(
    beta_hat = coefs, # GATE (point estimate)
    ghat_lower_point = tes_lci, # Pointwise Confidence Interval (lower bound)
    ghat_upper_point = tes_uci, # Pointwise Confidence Interval (upper bound)
    ghat_lower = tes_lcb, # Simultaneous Confidence Interval (lower bound)
    ghat_upper = tes_ucb, # Simultaneous Confidence Interval (upper bound)
    crit_val = crit_val
  )
}




