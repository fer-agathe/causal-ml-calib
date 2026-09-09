# ______________________________________________________________________________
# Model Evaluation & Diagnostic Metrics
# ______________________________________________________________________________

#' Expected Calibration Error (ECE)
#'
#' @description
#' Calculates the Expected Calibration Error (ECE) using uniform-width or quantile binning.
#'
#' @param d Binary indicator vector of observed labels/treatments.
#' @param scores Continuous vector of estimated predicted probabilities.
#' @param n_bins Integer specifying number of calibration bins. Default is `10`.
#' @param binning Character string indicating bin construction: `"width"` (equal length)
#'   or `"quantile"` (equal frequency). Default is `"width"`.
#'
#' @returns Numeric scalar value representing ECE.
#'
#' @md
compute_ece <- function(d,
                        scores,
                        n_bins = 10,
                        binning = "width") {
  n <- length(scores)

  if (binning == "width") {
    bin_cuts <- seq(0, 1, length.out = n_bins + 1)
  } else if (binning == "quantile") {
    bin_cuts <- quantile(
      scores, probs = seq(0, 1, length.out = n_bins + 1), names = FALSE
    )
    bin_cuts <- unique(bin_cuts)
  } else {
    stop("`binning` must be either 'width' or 'quantile'.")
  }

  # Assign bins
  bins <- cut(scores, breaks = bin_cuts, include.lowest = TRUE, labels = FALSE)

  # Drop empty bins
  keep <- !is.na(bins)
  bins <- bins[keep]
  d <- d[keep]
  scores <- scores[keep]

  # Aggregate statistics
  bin_counts <- as.numeric(table(bins))
  weights <- bin_counts / sum(bin_counts)

  acc <- tapply(d, bins, mean)
  if (any(is.na(acc))) acc[which(is.na(acc))] <- 0
  conf <- tapply(scores, bins, mean)
  if (any(is.na(conf))) conf[which(is.na(conf))] <- 0

  sum(weights * abs(acc - conf))
}


#' Kullback-Leibler (KL) Divergence
#'
#' @description
#' Calculates KL divergence between true underlying probability distribution $P$ and
#' estimated propensity score distribution $Q$.
#'
#' @param true_probas Numeric vector of ground-truth probabilities ($P$).
#' @param scores Numeric vector of predicted/calibrated probabilities ($Q$).
#' @param eps Small constant added for numerical stability to prevent $\log(0)$. Default is `1e-12`.
#'
#' @returns Numeric scalar value representing mean KL divergence.
#'
#' @md
kullback_leibler <- function(true_probas, scores, eps = 1e-12) {
  p <- pmin(pmax(true_probas, eps), 1 - eps)
  q <- pmin(pmax(scores, eps), 1 - eps)

  # Pointwise Bernoulli KL divergence: P log(P/Q) + (1-P) log((1-P)/(1-Q))
  kl_i <- p * log(p / q) + (1 - p) * log((1 - p) / (1 - q))
  mean(kl_i)
}
