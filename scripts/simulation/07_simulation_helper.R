# ______________________________________________________________________________
# Core functions and metrics computation for Monte Carlo simulation workflow
# ______________________________________________________________________________

library(tidyverse)
library(gmish)

#' Compute Performance Metrics for Propensity Score and Outcome Models
#'
#' @description
#' Evaluates error metrics and calibration performance for estimated nuisance
#' models within a single simulation iteration.
#'
#' @param y Numeric vector of observed outcomes.
#' @param d Numeric vector of binary treatment indicators (\eqn{D \in \{0, 1\}}).
#' @param m_true Numeric vector of ground-truth propensity scores.
#' @param m_pred Numeric vector of estimated (and potentially calibrated) propensity scores.
#' @param g0_pred Numeric vector of estimated potential outcomes under control (\eqn{\hat{g}_0(X)}).
#' @param g1_pred Numeric vector of estimated potential outcomes under treatment (\eqn{\hat{g}_1(X)}).
#' @param g_name Character string specifying the outcome learner used.
#' @param m_name Character string specifying the propensity score learner used.
#' @param alg_info Single-row data frame/tibble containing `algorithm_name` and `calib_name`.
#'
#' @return A single-row tibble containing RMSE for \eqn{g_0}, \eqn{g_1},
#' \eqn{m}, Brier score, ECE (width & quantile), Kullback-Leibler divergence,
#' and score bounds.
#'
#' @md
compute_simul_metrics <- function(y,
                                  d,
                                  m_true,
                                  m_pred,
                                  g0_pred,
                                  g1_pred,
                                  g_name,
                                  m_name,
                                  alg_info) {
  tibble(
    g_name         = g_name,
    m_name         = m_name,
    algorithm_name = alg_info$algorithm_name,
    calib_name     = alg_info$calib_name,
    # Metrics for outcome models
    rmse_g0        = sqrt(mean((y - g0_pred)^2)),
    rmse_g1        = sqrt(mean((y - g1_pred)^2)),
    # Metrics for propensity score models
    rmse_m         = sqrt(mean((m_true - m_pred)^2)),
    brier_m        = gmish::brier(pred = m_pred, obs = d),
    ece_width_m    = compute_ece(d = d, scores = m_pred, n_bins = 10, binning = "width"),
    ece_quantile_m = compute_ece(d = d, scores = m_pred, n_bins = 10, binning = "quantile"),
    kl_m           = kullback_leibler(true_probas = m_true, scores = m_pred),
    min_m          = min(m_pred),
    max_m          = max(m_pred)
  )
}


#' Compute Group Average Treatment Effects (GATE)
#'
#' @description
#' Constructs orthogonalized doubly robust signals and computes group average
#' treatment effects conditioned on a covariate feature.
#'
#' @param y Numeric vector of observed outcomes.
#' @param d Numeric vector of binary treatment indicators (\eqn{D \in \{0, 1\}}).
#' @param x_g Numeric vector of the conditioning covariate used for grouping.
#' @param g0_pred Numeric vector of estimated potential outcomes under control (\eqn{\hat{g}_0(X)}).
#' @param g1_pred Numeric vector of estimated potential outcomes under treatment (\eqn{\hat{g}_1(X)}).
#' @param m_pred Numeric vector of estimated propensity scores.
#' @param g_name Character string specifying the outcome learner used.
#' @param m_name Character string specifying the propensity score learner used.
#' @param alg_info Single-row data frame/tibble containing `algorithm_name` and `calib_name`.
#'
#' @return A single-row tibble containing estimated GATE parameters.
#'
#' @export
compute_simul_gate <- function(y,
                               d,
                               x_g,
                               g0_pred,
                               g1_pred,
                               m_pred,
                               g_name,
                               m_name,
                               alg_info) {
  # Construct orthogonal doubly robust signal
  robust_signal <- (y - g1_pred) * d / m_pred -
    (y - g0_pred) * (1 - d) / (1 - m_pred) +
    g1_pred - g0_pred

  gate_res <- group_average_treatment_effect(X = x_g, Y = robust_signal)

  as_tibble_row(gate_res$beta_hat) |>
    rename_with(~str_replace(.x, "Xraw", "gate"), .cols = everything()) |>
    mutate(
      g_name         = g_name,
      m_name         = m_name,
      algorithm_name = alg_info$algorithm_name,
      calib_name     = alg_info$calib_name,
      .before        = "gate1"
    )
}


#' Resolve Clipping Threshold and Calibration Specification
#'
#' @description
#' Helper function to resolve dynamic clipping bounds (derived from Platt or Beta
#' calibration ranges in preceding loops) or default static thresholds.
#'
#' @param calib_row Single-row data frame specifying `algorithm_name` and `calib_name`.
#' @param g_name Character string specifying outcome model learner name.
#' @param m_name Character string specifying propensity score learner name.
#' @param range_scores Tibble recording score bounds from previous calibration fits.
#' @param default_clipping Baseline numeric threshold or 2-element numeric vector.
#'
#' @return A list with elements `clip_threshold` and `calib_type`.
#'
#' @export
resolve_clipping_and_calib <- function(calib_row,
                                       g_name,
                                       m_name,
                                       range_scores,
                                       default_clipping) {
  if (calib_row$calib_name %in% c("uncalib-platt-clip", "uncalib-beta-clip")) {
    target_calib <- ifelse(str_detect(calib_row$calib_name, "platt"), "platt", "beta")

    range_calib <- range_scores |>
      filter(
        g_name == !!g_name,
        m_name == !!m_name,
        algorithm_name == calib_row$algorithm_name,
        calib_name == !!target_calib
      )

    list(
      clip_threshold = c(range_calib$min, range_calib$max),
      calib_type     = "uncalibrated"
    )
  } else {
    list(
      clip_threshold = default_clipping,
      calib_type     = calib_row$calib_name
    )
  }
}


#' Single Monte Carlo Simulation Iteration
#'
#' @description
#' Runs a single simulation replicate for a specified DGP. Fits multiple outcome
#' and propensity score models, applies calibration routines across $K$-folds,
#' and evaluates ATE, GATE, and model quality metrics.
#'
#' @param dgp_name Character string identifying the DGP function (e.g.,
#'   `"belloni"`).
#' @param dgp_params List of parameters passed to the DGP data generation function.
#' @param clipping_threshold Default threshold or range for propensity score clipping.
#'   Default is `1e-12`.
#' @param n_folds Integer specifying the number of folds for cross-fitting.
#'   Default is `5`.
#' @param equilibrium Logical indicating if outcome models account for
#'   equilibrium conditions. Default is `TRUE`.
#' @param seed Optional integer to set pseudo-random state for reproducible replication.
#'
#' @details
#' Algorithms supported include:
#'   *`alg-1-uncalibrated`: Standard DML without propensity score calibration.
#'   *`alg-3-cross-fitted-calib`: DML incorporating cross-fitted propensity
#'     score calibration (Isotonic, Platt, Beta, or custom clip).
#'
#' @return A named list containing:
#' * `ate`: Tibble of estimated ATE values across models and algorithms.
#' * `gate`: Tibble of estimated GATE parameters across models and algorithms.
#' * `metrics`:  Tibble of error and calibration metrics for nuisance models.
#' * `dgp_params`: Input DGP parameter list.
#' * `range_scores`: Tibble recording minimum and maximum predicted propensity scores.
#' * `clipping_threshold`: Input clipping threshold configuration.
#'
#' @md
simul <- function(dgp_name,
                  dgp_params,
                  clipping_threshold = 1e-12,
                  n_folds = 5,
                  equilibrium = TRUE,
                  seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # 1. Data Generation & Formatting ----
  if (dgp_name == "belloni") {
    data_raw <- gen_data_belloni(
      n   = dgp_params$n,
      p   = dgp_params$p,
      Rd2 = dgp_params$Rd2,
      Ry2 = dgp_params$Ry2,
      rho = dgp_params$rho,
      seed = seed
    )
    tb_data   <- format_data_belloni(data_raw)
    ate_true  <- data_raw$ate_true
    gate_true <- data_raw$gate_true
  } else {
    stop(sprintf("Unsupported DGP specification: '%s'", dgp_name))
  }

  # Method Configurations
  calib_methods <- tribble(
    ~algorithm_name,            ~calib_name,
    'alg-1-uncalibrated',       'uncalibrated',
    'alg-3-cross-fitted-calib', 'isotonic',
    'alg-3-cross-fitted-calib', 'platt',
    'alg-3-cross-fitted-calib', 'uncalib-platt-clip',
    'alg-3-cross-fitted-calib', 'beta',
    'alg-3-cross-fitted-calib', 'uncalib-beta-clip',
    'alg-3-cross-fitted-calib', 'joint-platt'
  )

  g_names <- c("ols", "lasso", "rf", "lgbm")
  m_names <- c("logit", "rf_classif", "rf_reg", "lgbm")

  # Fold Partitioning
  folds     <- sample(rep(1:n_folds, length.out = nrow(tb_data)))
  ind_folds <- map(seq_len(n_folds), ~which(folds == .x)) |> list_c()

  # 2. Fit Outcome & Propensity Score Learners ----
  models_g <- map(
    g_names,
    ~estim_g(
      data = tb_data, learner_g_name = .x, folds = folds, equilibrium = equilibrium
    )) |> set_names(g_names)
  models_m <- map(
    m_names, ~estim_m(data = tb_data, learner_m_name = .x, folds = folds)
    ) |> set_names(m_names)

  # Result Storage Accumulators
  ate_list          <- list()
  gate_list         <- list()
  metrics_list      <- NULL
  range_scores <- NULL
  # calib_curves_list <- NULL

  # 3. Main Evaluation Loops ----
  for (g_name in g_names) {
    for (m_name in m_names) {
      for (i_algo in seq_len(nrow(calib_methods))) {

        alg_info <- calib_methods[i_algo, ]

        # Determine current clipping threshold and calibration settings
        calib_cfg <- resolve_clipping_and_calib(
          calib_row        = alg_info,
          g_name           = g_name,
          m_name           = m_name,
          range_scores     = range_scores,
          default_clipping = clipping_threshold
        )

        # Estimate Nuisances
        nuisance_curr <- nuisance_estim(
          g_name             = g_name,
          m_name             = m_name,
          algorithm_name     = alg_info$algorithm_name,
          calib_name         = calib_cfg$calib_type,
          data               = tb_data,
          folds              = folds,
          n_folds            = n_folds,
          models_g           = models_g,
          models_m           = models_m,
          clipping_threshold = calib_cfg$clip_threshold
        )

        # Extract predictions across folds
        g0_pred <- map(nuisance_curr, "g0_pred_k") |> list_c()
        g1_pred <- map(nuisance_curr, "g1_pred_k") |> list_c()
        geqb_pred <- map(nuisance_curr, "geqb_pred_k") |> list_c()
        m_pred  <- map(nuisance_curr, "m_pred_k") |> list_c()

        # Observed targets & true values
        y      <- tb_data$y[ind_folds]
        d      <- tb_data$d[ind_folds]
        m_true <- data_raw$m_0[ind_folds]

        # ATE Computation
        ate_curr <- compute_ate_folds(pred_folds = nuisance_curr, n_folds = n_folds) |>
          mutate(
            g_name         = g_name,
            m_name         = m_name,
            algorithm_name = alg_info$algorithm_name,
            calib_name     = alg_info$calib_name
          )
        ate_list[[length(ate_list) + 1]] <- ate_curr

        # GATE Computation
        gate_curr <- compute_simul_gate(
          y = y, d = d, x_g = tb_data$X1[ind_folds],
          g0_pred = g0_pred, g1_pred = g1_pred, m_pred = m_pred,
          g_name = g_name, m_name = m_name, alg_info = alg_info
        )
        gate_list[[length(gate_list) + 1]] <- gate_curr

        # Metrics Computation
        metrics_curr <- compute_simul_metrics(
          y = y, d = d, m_true = m_true, m_pred = m_pred,
          g0_pred = g0_pred, g1_pred = g1_pred,
          g_name = g_name, m_name = m_name, alg_info = alg_info
        )
        metrics_list[[length(metrics_list) + 1]] <- metrics_curr

        # Track predicted range scores
        range_scores_curr <- tibble(
          g_name         = g_name,
          m_name         = m_name,
          algorithm_name = alg_info$algorithm_name,
          calib_name     = alg_info$calib_name,
          min            = metrics_curr$min_m,
          max            = metrics_curr$max_m
        )
        range_scores <- range_scores |>
          bind_rows(range_scores_curr)

        # # Calibration curve
        # tb_current <- tibble(
        #   id = ind_folds,
        #   y = y,
        #   d = d,
        #   m_true = m_true,
        #   m_pred = m_pred
        # ) |>
        #   mutate(
        #     g_name = g_name,
        #     m_name = m_name,
        #     algorithm_name = alg_info$algorithm_name,
        #     calib_name = alg_info$calib_name
        #   )
        # calib_pred <- calibration_curve(data = tb_current, d_name = "d", p = m_pred)
        # calib_curves_curr <- list(
        #   g_name = g_name,
        #   m_name = m_name,
        #   algorithm_name = alg_info$algorithm_name,
        #   calib_name = alg_info$calib_name,
        #   calib_pred = calib_pred
        # )
        # calib_curves_list[[length(calib_curves_list) + 1]] <- calib_curves_curr
      }
    }
  }

  # 4. Final Output Construction ----
  list(
    ate = bind_rows(ate_list) |>
      mutate(seed = seed, ate_true = ate_true),
    gate = bind_rows(gate_list) |>
      mutate(
        seed       = seed,
        gate_true1 = gate_true[1],
        gate_true2 = gate_true[2],
        gate_true3 = gate_true[3],
        gate_true4 = gate_true[4],
        gate_true5 = gate_true[5]
      ),
    metrics            = bind_rows(metrics_list) |> mutate(seed = seed),
    dgp_params         = dgp_params,
    range_scores       = range_scores,
    # calib_curves       = calib_curves_list,
    clipping_threshold = clipping_threshold
  )
}
