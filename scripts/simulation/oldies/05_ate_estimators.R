# ______________________________________________________________________________
# Average Treatment Effect (ATE) Estimators
# ______________________________________________________________________________

#' Compute AIPW / Doubly Robust ATE Estimator
#'
#' Computes the Augmented Inverse Probability Weighting (AIPW) doubly robust estimator
#' for the Average Treatment Effect (ATE), alongside its asymptotic standard error
#' and confidence interval.
#'
#' @details
#' The doubly robust ATE score for unit \eqn{i} is given by:
#' \deqn{\psi_i = \hat{g}_1(X_i) - \hat{g}_0(X_i) + \frac{D_i (Y_i - \hat{g}_1(X_i))}{\hat{m}(X_i)} - \frac{(1 - D_i)(Y_i - \hat{g}_0(X_i))}{1 - \hat{m}(X_i)}}
#' The point estimate is the sample average \eqn{\hat{\tau}_{\text{AIPW}} = \frac{1}{n} \sum_{i=1}^n \psi_i}.
#'
#' @param y Numeric vector of observed outcomes \eqn{Y}.
#' @param d Integer or numeric vector of binary treatment indicators \eqn{D \in \{0, 1\}}.
#' @param g0_hat Numeric vector of predicted baseline outcome values \eqn{\hat{g}_0(X) = \mathbb{E}[Y \mid D=0, X]}.
#' @param g1_hat Numeric vector of predicted treated outcome values \eqn{\hat{g}_1(X) = \mathbb{E}[Y \mid D=1, X]}.
#' @param m_hat Numeric vector of predicted or calibrated propensity scores \eqn{\hat{m}(X) = \Pr(D=1 \mid X)}.
#' @param alpha Significance level for the nominal \eqn{(1 - \alpha) \times 100\%} confidence interval.
#'   Default is `0.05`.
#'
#' @returns A tibble with a single row containing:
#'   * `estimate`: Point estimate of the ATE \eqn{\hat{\tau}_{\text{AIPW}}}.
#'   * `se`: Asymptotic standard error derived from the influence function.
#'   * `ci_lower`: Lower bound of the nominal \eqn{(1 - \alpha)} confidence interval.
#'   * `ci_upper`: Upper bound of the nominal \eqn{(1 - \alpha)} confidence interval.
#'
#' @md
ate_aipw <- function(y, d, g0_hat, g1_hat, m_hat, alpha = 0.05) {
  n <- length(y)

  psi <- (g1_hat - g0_hat) +
    (d * (y - g1_hat) / m_hat) -
    ((1 - d) * (y - g0_hat) / (1 - m_hat))

  estimate <- mean(psi)
  se <- stats::sd(psi) / sqrt(n)
  z_crit <- stats::qnorm(1 - alpha / 2)

  tibble::tibble(
    estimate = estimate,
    se = se,
    ci_lower = estimate - z_crit * se,
    ci_upper = estimate + z_crit * se
  )
}


#' Compute Inverse Probability Weighting (IPW) ATE Estimator
#'
#' Computes the standard Inverse Probability Weighting (IPW) estimator for the ATE
#' without regression adjustment.
#'
#' @details
#' The IPW score for observation \eqn{i} is defined as:
#' \deqn{\psi_i^{\text{IPW}} = \frac{D_i Y_i}{\hat{m}(X_i)} - \frac{(1 - D_i) Y_i}{1 - \hat{m}(X_i)}}
#'
#' @param y Numeric vector of observed outcomes \eqn{Y}.
#' @param d Integer or numeric vector of binary treatment indicators \eqn{D \in \{0, 1\}}.
#' @param m_hat Numeric vector of predicted propensity scores \eqn{\hat{m}(X) = \Pr(D=1 \mid X)}.
#' @param alpha Significance level for the confidence interval. Default is `0.05`.
#'
#' @returns A tibble containing `estimate`, `se`, `ci_lower`, and `ci_upper`.
#'
#' @md
ate_ipw <- function(y,
                    d,
                    m_hat,
                    alpha = 0.05) {
  n <- length(y)

  psi <- (d * y / m_hat) - ((1 - d) * y / (1 - m_hat))

  estimate <- mean(psi)
  se <- stats::sd(psi) / sqrt(n)
  z_crit <- stats::qnorm(1 - alpha / 2)

  tibble::tibble(
    estimate = estimate,
    se = se,
    ci_lower = estimate - z_crit * se,
    ci_upper = estimate + z_crit * se
  )
}


#' Compute Outcome Regression (OR) ATE Estimator
#'
#' Computes the simple Outcome Regression (G-computation) estimator for the ATE.
#'
#' @details
#' The outcome regression estimator is given by the average difference in
#' conditional potential outcome predictions:
#' \deqn{\hat{\tau}_{\text{OR}} = \frac{1}{n} \sum_{i=1}^n \left( \hat{g}_1(X_i) - \hat{g}_0(X_i) \right)}
#'
#' @param g0_hat Numeric vector of predicted baseline outcomes \eqn{\hat{g}_0(X)}.
#' @param g1_hat Numeric vector of predicted treated outcomes \eqn{\hat{g}_1(X)}.
#' @param alpha Significance level for the confidence interval. Default is `0.05`.
#'
#' @returns A tibble containing `estimate`, `se`, `ci_lower`, and `ci_upper`.
#'
#' @md
ate_reg <- function(g0_hat,
                    g1_hat,
                    alpha = 0.05) {
  n <- length(g0_hat)

  psi <- g1_hat - g0_hat

  estimate <- mean(psi)
  se <- stats::sd(psi) / sqrt(n)
  z_crit <- stats::qnorm(1 - alpha / 2)

  tibble::tibble(
    estimate = estimate,
    se = se,
    ci_lower = estimate - z_crit * se,
    ci_upper = estimate + z_crit * se
  )
}


#' Compare All ATE Estimators
#'
#' Wrapper function to run AIPW, IPW, and Outcome Regression ATE estimators
#' simultaneously on a single set of cross-fitted nuisance predictions.
#'
#' @param y Numeric vector of observed outcomes \eqn{Y}.
#' @param d Integer or numeric vector of binary treatment indicators \eqn{D}.
#' @param nuisance_preds A tibble or list containing predicted nuisance
#'   columns `g0_hat`, `g1_hat`, and `m_hat`, as returned by `\link{predict_nuisance_calibrated}`.
#' @param alpha Significance level for confidence intervals. Default is `0.05`.
#'
#' @returns A tibble with one row per estimator (`"AIPW"`, `"IPW"`, `"Regression"`):
#'  * `estimator`: Character string identifying the estimation method.
#'  * `estimate`: Point estimate of the ATE.
#'  * `se`: Standard error of the estimate.
#'  * `ci_lower`: Lower bound of the \eqn{(1 - \alpha)} confidence interval.
#'  * `ci_upper`: Upper bound of the \eqn{(1 - \alpha)} confidence interval.
#'
#' @md
estimate_all_ate <- function(y,
                             d,
                             nuisance_preds,
                             alpha = 0.05) {
  res_aipw <- ate_aipw(
    y = y,
    d = d,
    g0_hat = nuisance_preds$g0_hat,
    g1_hat = nuisance_preds$g1_hat,
    m_hat = nuisance_preds$m_hat,
    alpha = alpha
  ) |>
    dplyr::mutate(estimator = "AIPW")

  res_ipw <- ate_ipw(
    y = y,
    d = d,
    m_hat = nuisance_preds$m_hat,
    alpha = alpha
  ) |>
    dplyr::mutate(estimator = "IPW")

  res_reg <- ate_reg(
    g0_hat = nuisance_preds$g0_hat,
    g1_hat = nuisance_preds$g1_hat,
    alpha = alpha
  ) |>
    dplyr::mutate(estimator = "Regression")

  dplyr::bind_rows(res_aipw, res_ipw, res_reg) |>
    dplyr::select("estimator", "estimate", "se", "ci_lower", "ci_upper")
}
