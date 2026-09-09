# ______________________________________________________________________________
# Propensity Score Calibration, Nuisance Orchestration & Causal Inference
# ______________________________________________________________________________

library(tidyverse)
library(sandwich)
library(philentropy)
library(gmish)

#' Calibrate Predicted Probabilities / Propensity Scores
#'
#' @description
#' Fits a specified calibration model on validation/training predictions and
#' outputs calibrated propensity scores bounded within a clipping interval.
#'
#' @param d Numeric vector of binary treatment status (\eqn{D \in \{0, 1\}}).
#' @param scores Numeric vector of uncalibrated predicted probabilities
#'   (\eqn{\hat{m}(X) \in (0, 1)}).
#' @param new_scores Numeric vector of uncalibrated predicted probabilities on
#'   test data to be calibrated.
#' @param g0 Numeric vector of predicted outcomes in group 0 (\eqn{\hat{g}_0(X)}).
#'   Default to `NULL`, only used if `method="joint-platt`.
#' @param g0 Numeric vector of predicted outcomes in group 1 (\eqn{\hat{g}_0(X)}).
#'   Default to `NULL`, only used if `method="joint-platt`.
#' @param new_g0 Numeric vector of predicted outcomes in group 0 on test data.
#'   Default to `NULL`, only used if `method="joint-platt`.
#' @param new_g0 Numeric vector of predicted outcomes in group 1 on test data.
#'   Default to `NULL`, only used if `method="joint-platt`.
#' @param method Character string specifying calibration method. Options:
#'   `"uncalibrated"`, `"isotonic"`, `"platt"`, `"beta"`, or `"joint-platt"`.
#' @param clipping_threshold Numeric scalar or length-2 numeric vector
#'   specifying lower and upper clipping bounds. Default is `1e-12`
#'   (clipped to `[1e-12, 1 - 1e-12]`).
#'
#' @returns A numeric vector of calibrated predicted probabilities for `new_scores`.
#'
#' @md
calibrate <- function(d,
                      scores,
                      new_scores,
                      g0 = NULL,
                      g1 = NULL,
                      new_g0 = NULL,
                      new_g1 = NULL,
                      method = "uncalibrated",
                      clipping_threshold = 1e-12) {

  # Resolve upper/lower clipping boundaries
  if (length(clipping_threshold) == 1) {
    clip_lower <- clipping_threshold
    clip_upper <- 1 - clipping_threshold
  } else if (length(clipping_threshold) == 2) {
    clip_lower <- clipping_threshold[1]
    clip_upper <- clipping_threshold[2]
  } else {
    stop("`clipping_threshold` must be a numeric scalar or a length-2 numeric vector.")
  }

  if (method == "uncalibrated") {
    calib_pred <- new_scores

  } else {
    eps        <- 1e-12
    scores_safe <- pmin(pmax(scores, eps), 1 - eps)
    new_scores_safe   <- pmin(pmax(new_scores, eps), 1 - eps)

    if (method == "isotonic") {
      iso_fit    <- isoreg(x = scores_safe, y = d)
      step_fn <- as.stepfun(iso_fit)
      calib_pred <- step_fn(new_scores_safe)

    } else if (method == "platt") {
      platt_fit  <- glm(d ~ scores_safe, family = binomial(link = "logit"))
      calib_pred <- predict(
        platt_fit, newdata = data.frame(new_scores_safe = new_scores_safe),
        type = "response"
      )

    } else if (method == "beta") {
      log_p       <- log(scores_safe)
      log_1p      <- log(1 - scores_safe)
      new_log_p   <- log(new_scores_safe)
      new_log_1p  <- log(1 - new_scores_safe)

      beta_fit   <- glm(d ~ log_p + log_1p, family = binomial(link = "logit"))
      calib_pred <- predict(
        beta_fit, newdata = data.frame(log_p = new_log_p, log_1p = new_log_1p),
        type = "response"
      )

    } else if (method == "joint-platt") {

      feat_mat <- cbind(scores_safe, g0, g1)
      new_feat_mat <- cbind(new_scores_safe, g0, g1)

      feat_mat <- scale(feat_mat) # Standardize columns
      new_feat_mat <- scale(new_feat_mat)

      joint_platt_fit <- glm(
        d ~ ., data = as.data.frame(feat_mat), family = binomial(link = "logit")
      )
      calib_pred <- predict(
        joint_platt_fit,
        newdata = as.data.frame(new_feat_mat), type = "response")

    } else {
      stop(sprintf("Unsupported calibration method: '%s'", method))
    }
  }


  # Enforce probability clipping boundaries
  pmin(pmax(calib_pred, clip_lower), clip_upper)
}


#' Estimate Nuisance Parameters Across Folds (Algorithm Router)
#'
#' @description
#' High-level wrapper function to estimate cross-fitted outcome regressions and
#' calibrated propensity scores for a specific algorithm architecture.
#'
#' @param g_name Character string specifying outcome model learner name.
#' @param m_name Character string specifying propensity score learner name.
#' @param algorithm_name Character string identifying algorithm structure:
#'   `"alg-1-uncalibrated"` or `"alg-3-cross-fitted-calib"`.
#' @param calib_name Character string specifying the calibration method passed
#'   to \code{\link{calibrate}}.
#' @param data Prepared dataset (tibble/data.frame) containing variables `y`,
#'  `d`, and feature matrix columns.
#' @param folds Vector of fold assignment integers for each observation.
#' @param n_folds Integer specifying number of cross-fitting folds.
#' @param models_g Pre-fitted outcome models list produced by \code{estim_g()}.
#' @param models_m Pre-fitted propensity score models list produced by \code{estim_m()}.
#' @param clipping_threshold Lower/upper clipping thresholds passed to \code{\link{calibrate}}.
#'
#' @returns A list of length `n_folds`, where each element contains out-of-fold predicted
#'   vectors: `g0_pred_k`, `g1_pred_k`, and `m_pred_k`.
#'
#' @md
nuisance_estim <- function(g_name,
                           m_name,
                           algorithm_name,
                           calib_name,
                           data,
                           folds,
                           n_folds,
                           models_g,
                           models_m,
                           clipping_threshold = 1e-12) {

  if (algorithm_name == "alg-1-uncalibrated") {
    alg_1_uncalibrated(
      g_name = g_name, m_name = m_name, data = data, folds = folds,
      n_folds = n_folds, models_g = models_g, models_m = models_m,
      clipping_threshold = clipping_threshold
    )
  } else if (algorithm_name == "alg-3-cross-fitted-calib") {
    alg_3_cross_fitted_calib(
      g_name = g_name, m_name = m_name, calib_name = calib_name,
      data = data, folds = folds, n_folds = n_folds,
      models_g = models_g, models_m = models_m,
      clipping_threshold = clipping_threshold
    )
  } else {
    stop(sprintf("Unknown algorithm structure: '%s'", algorithm_name))
  }
}


#' Standard Uncalibrated Double Machine Learning Nuisance Estimation (Algorithm 1)
#'
#' @description
#' Computes out-of-fold outcome and propensity predictions using standard K-fold
#' cross-fitting without additional calibration steps.
#'
#' @inheritParams nuisance_estim
#'
#' @returns A list of length `n_folds` containing out-of-fold predictions.
#'
#' @md
alg_1_uncalibrated <- function(g_name,
                               m_name,
                               data,
                               folds,
                               n_folds,
                               models_g,
                               models_m,
                               clipping_threshold = 1e-12) {

  map(seq_len(n_folds), function(k) {
    test_idx <- which(folds == k)
    data_k   <- data[test_idx, ]
    d_k <- data_k$d
    y_k <- data_k$y
    data_k <- data_k |> dplyr::select(-d, -y)

    # Predict outcome regression g(X)
    predict_g_fn <- get(sprintf("predict_g_%s", g_name))
    g_model_k    <- models_g[[g_name]][[k]] # previously trained outcome learner model (on I_k^c)
    g0_preds     <- predict_g_fn(g_model_k$fitted_g_0, data_k)
    g1_preds     <- predict_g_fn(g_model_k$fitted_g_1, data_k)
    if (!is.null(g_model_k$fitted_g_eqb)) {
      geqb_preds <- predict_g_fn(g_model_k$fitted_g_eqb, data_k)
    } else geqb_pred_k <- NULL

    # Predict propensity score m(X)
    predict_m_fn <- get(sprintf("predict_m_%s", m_name))
    m_raw        <- predict_m_fn(models_m[[m_name]][[k]], data_k)

    # Clip raw propensity predictions
    m_clipped    <- calibrate(
      d                  = d_k,
      scores             = m_raw,
      new_scores         = m_raw,
      g0                 = g0_preds,
      g1                 = g1_preds,
      new_g0             = g0_preds,
      new_g1             = g1_preds,
      method             = "uncalibrated",
      clipping_threshold = clipping_threshold
    )

    list(
      k           = k,
      g0_pred_k   = g0_preds,
      g1_pred_k   = g1_preds,
      geqb_pred_k = geqb_preds,
      m_pred_k    = m_clipped,
      y_k         = y_k,
      d_k         = d_k
    )
  })
}


#' Cross-Fitted Calibrated Double Machine Learning Nuisance Estimation (Algorithm 3)
#'
#' @description
#' Performs cross-fitted calibration where primary nuisance predictions from
#' fold k are calibrated on fold k' before generating final out-of-fold
#' predictions for target fold k.
#'
#' @inheritParams nuisance_estim
#'
#' @returns A list of length `n_folds` containing out-of-fold predictions.
#'
#' @md
alg_3_cross_fitted_calib <- function(g_name,
                                     m_name,
                                     calib_name,
                                     data,
                                     folds,
                                     n_folds,
                                     models_g,
                                     models_m,
                                     clipping_threshold = 1e-12) {

  predict_g_fn <- get(sprintf("predict_g_%s", g_name))
  predict_m_fn <- get(sprintf("predict_m_%s", m_name))

  map(seq_len(n_folds), function(k) {
    test_idx <- which(folds == k)
    data_k   <- data[test_idx, ]
    d_k <- data_k$d
    y_k <- data_k$y
    data_k <- data_k |> dplyr::select(-d, -y)

    # Primary predictions on test fold k
    g_model_k <- models_g[[g_name]][[k]] # previously trained outcome learner model (on I_k^c)
    g0_preds <- predict_g_fn(g_model_k$fitted_g_0, data_k)
    g1_preds <- predict_g_fn(g_model_k$fitted_g_1, data_k)
    if (!is.null(g_model_k$fitted_g_eqb)) {
      geqb_preds <- predict_g_fn(g_model_k$fitted_g_eqb, data_k)
    } else geqb_preds <- NULL

    m_raw_k <- predict_m_fn(models_m[[m_name]][[k]], data_k)

    # Clipping
    m_clipped_k <- pmin(
      pmax(m_raw_k, 0.0 + clipping_threshold),
      1.0 - clipping_threshold
    )

    # Calibration based on data from the left-out sample
    m_calibrated <- calibrate(
      d                  = d_k,
      scores             = m_clipped_k,
      new_scores         = m_clipped_k,
      g0                 = g0_preds,
      g1                 = g1_preds,
      new_g0             = g0_preds,
      new_g1             = g1_preds,
      method             = calib_name,
      clipping_threshold = clipping_threshold
    )

    list(
      k           = k,
      g0_pred_k   = g0_preds,
      g1_pred_k   = g1_preds,
      geqb_pred_k = geqb_preds,
      m_pred_k    = m_calibrated,
      y_k         = y_k,
      d_k         = d_k
    )
  })
}



calibration_curve <- function(data, d_name, p){
  D <- data |> pull(!!d_name)
  cls <- cut(p, breaks = unique(c(0, quantile(p, (1:9) / 10), 1)))
  agg_y <- aggregate(D, by = list(cls), mean)
  agg_p <- aggregate(p, by = list(cls), mean)
  list(
    level = agg_y$Group.1,
    score = agg_p$x,
    y     = agg_y$x
  )
}
