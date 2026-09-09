# ______________________________________________________________________________
# Monte Carlo Simulation
# ______________________________________________________________________________

# Source modular components
source("scripts/simulation/01_dgp.R")
source("scripts/simulation/02_learners_g.R")
source("scripts/simulation/03_learners_m.R")
source("scripts/simulation/04_calibration.R")
source("scripts/simulation/05_metrics.R")
source("scripts/simulation/06_causal-estimators.R")
source("scripts/simulation/07_simulation_helper.R")

library(tidyverse)
library(pbapply)
library(parallel)

# 1. Parameter Grid Setup ----


# ~3h for each simulation with 20 variables
# ~18h for each simulation with 200 variables
grid <- expand_grid(
  n    = c(2000),
  p    = c(20, 200),
  Rd2  = c(0.2, 0.5, 0.8),
  clip = 1e-12 # 1e-12 corresponds to virtually no clipping
)

output_dir <- "output/simul-klassen/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 2. Parallel Setup ----

nb_simul <- 1000
ncl <- max(1, detectCores() - 1)
cl  <- makeCluster(ncl)

# Guarantee cluster shutdown even if an error occurs mid-script
on.exit(stopCluster(cl), add = TRUE)

# Load dependencies across all worker nodes
clusterEvalQ(cl, {
  library(tidyverse)
  library(hdm)
  library(randomForest)
  library(xgboost)
  library(sandwich)
  library(philentropy)
}) |> invisible()

# Export custom helper and estimation functions
functions_to_export <- c(
  "simul",
  # Data Generation
  "gen_data_belloni", "format_data_belloni",
  # Learners g
  "estim_g", "estim_g_ols", "predict_g_ols",
  "estim_g_lasso", "predict_g_lasso",
  "estim_g_rf", "predict_g_rf",
  "estim_g_lgbm", "predict_g_lgbm",
  # Learners m
  "estim_m",
  "estim_m_logit", "predict_m_logit",
  "estim_m_rf_classif", "predict_m_rf_classif",
  "estim_m_rf_classif", "predict_m_rf_classif",
  "estim_m_rf_reg", "predict_m_rf_reg",
  "estim_m_lgbm", "predict_m_lgbm",
  # Calibration
  "calibrate", "nuisance_estim", "alg_1_uncalibrated",
  "alg_3_cross_fitted_calib", "calibration_curve",
  # Metrics
  "compute_ece", "kullback_leibler",
  # Causal Estimators
  "compute_ate_reg", "compute_ate_ipw", "compute_ate_aipw",
  "compute_ate_ipw_normalized", "compute_ate_aipw_normalized",
  "compute_ate_eqb_reg", "compute_ate_eqb_ipw", "compute_ate_eqb_aipw",
  "compute_ate_folds",
  "qtmax", "group_average_treatment_effect",
  # Simulation Helpers
  "compute_simul_metrics", "compute_simul_gate", "resolve_clipping_and_calib"
)

clusterExport(cl, functions_to_export)

# 3. Main Simulation Loop ----

cat(sprintf("Starting Monte Carlo simulations (%d combinations across %d cores)...\n\n",
            nrow(grid), ncl))

for (i_grid in 1:nrow(grid)) {
  cat(
    sprintf(
      "## Simulation %d/%d (n=%d, p=%d, Rd2=%.1f, clip=%g)\n",
      i_grid, nrow(grid), grid$n[i_grid], grid$p[i_grid],
      grid$Rd2[i_grid], grid$clip[i_grid])
  )

  dgp_params <- list(
    n   = grid$n[i_grid],
    p   = grid$p[i_grid],
    Rd2 = grid$Rd2[i_grid],
    Ry2 = 0.5,
    rho = 0.5
  )

  clipping_threshold <- grid$clip[i_grid]

  # Export current parameter state to worker nodes
  clusterExport(cl, c("dgp_params", "clipping_threshold"), envir = environment())

  # Run replicates in parallel
  res_tmp <- pblapply(1:nb_simul, FUN = function(seed) {
    simul(
      dgp_name           = "belloni",
      dgp_params         = dgp_params,
      clipping_threshold = clipping_threshold,
      n_folds            = 5,
      seed               = seed
    )
  }, cl = cl)

  # Save individual grid output iteration
  save(res_tmp, file = str_c(output_dir, i_grid, ".rda"))
}

cat("\nSimulations completed successfully.\n")

# 4. Process and Aggregate Results ------------------------------------------

cat("Processing and assembling results structure...\n")

res <- list()

for (i_grid in 1:nrow(grid)) {
  cat(
    sprintf(
      "Load results from simulation (%d/%d)...\n\n",
      i_grid, nrow(grid)
    )
  )
  file_path <- str_c(output_dir, i_grid, ".rda")

  if (!file.exists(file_path)) next

  load(file_path) # Loads 'res_tmp'

  for (j in seq_along(res_tmp)) {
    # Extract DGP parameters
    cur_dgp  <- res_tmp[[j]]$dgp_params
    cur_clip <- res_tmp[[j]]$clipping_threshold

    # Annotate ATE table
    res_tmp[[j]]$ate <- res_tmp[[j]]$ate |>
      mutate(
        n                  = cur_dgp$n,
        p                  = cur_dgp$p,
        Rd2                = cur_dgp$Rd2,
        Ry2                = cur_dgp$Ry2,
        rho                = cur_dgp$rho,
        clipping_threshold = cur_clip
      )

    # Annotate GATE table
    res_tmp[[j]]$gate <- res_tmp[[j]]$gate |>
      mutate(
        n                  = cur_dgp$n,
        p                  = cur_dgp$p,
        Rd2                = cur_dgp$Rd2,
        Ry2                = cur_dgp$Ry2,
        rho                = cur_dgp$rho,
        clipping_threshold = cur_clip
      )

    # Annotate Metrics table
    res_tmp[[j]]$metrics <- res_tmp[[j]]$metrics |>
      mutate(
        n                  = cur_dgp$n,
        p                  = cur_dgp$p,
        Rd2                = cur_dgp$Rd2,
        Ry2                = cur_dgp$Ry2,
        rho                = cur_dgp$rho,
        clipping_threshold = cur_clip
      )
  }

  res <- c(res, res_tmp)
}

# Combine into master dataframes
res_ate     <- map(res, "ate") |> list_rbind()
res_gate    <- map(res, "gate") |> list_rbind()
res_metrics <- map(res, "metrics") |> list_rbind()

# 5. Summarize Metrics ----

metrics_ate <- res_ate |>
  mutate(
    delta_ate_reg  = ate_reg - ate_true,
    delta_ate_ipw  = ate_ipw - ate_true,
    delta_ate_aipw = ate_aipw - ate_true
  ) |>
  group_by(g_name, m_name, algorithm_name, calib_name, n, p, Rd2, Ry2, rho, clipping_threshold) |>
  summarise(
    mae_ate_reg   = mean(abs(delta_ate_reg)),
    mae_ate_ipw   = mean(abs(delta_ate_ipw)),
    mae_ate_aipw  = mean(abs(delta_ate_aipw)),
    rmse_ate_reg  = sqrt(mean(delta_ate_reg^2)),
    rmse_ate_ipw  = sqrt(mean(delta_ate_ipw^2)),
    rmse_ate_aipw = sqrt(mean(delta_ate_aipw^2)),
    sd_ate_reg    = sd(ate_reg),
    sd_ate_ipw    = sd(ate_ipw),
    sd_ate_aipw   = sd(ate_aipw),
    .groups       = "drop"
  )


