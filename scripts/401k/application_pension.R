# _______________________________________________________________
# Empirical illustration: the 401(k) data
# _______________________________________________________________
# 1. Semi-synthetic experiment: known ATE
# __
# For each replication, draw a subsample of n=2000 households
# learners receive 9 real covariates and 191 noise proxies, leading to p=200.
#
# Output: output/401k/results_mc_pension.csv
# rep, model, calib, estimator, estimate, ece, r2g0, ATE_true
#
# Duration: ~10 minutes for 1,000 replications
#
# 2. Actual data
# __
# Using the whole e401 dataset:
# net_tfa, AIPW (uncal/Hájek/Beta) + IPW/IPW-N, SE
#
# Output: output/401k/estimates_real_pension.csv
# model, calib, estimator, est, se
#
# Duration: ~30 seconds
# _______________________________________________________________

output_dir <- "output/401k/"
fig_dir <- "output/figs/"
tbl_dir <- "output/tbl"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(fig_dir)) dir.create(output_dir, recursive = TRUE)

# Turn to TRUE to export graphs in PDF
export_graphs <- TRUE

# 1 Settings----

# Themes for ggplot2 graphs
source("scripts/theme.R")

library(hdm)
library(ranger)
library(xgboost)
library(tidyverse)

set.seed(20260712)
NREP  <- 1000  # Monte Carlo replications
NSUB  <- 2000  # Size of sub-samples
P_ADD <- 191   # No. added noisy proxies (p = 9 + P_ADD)
KF_MC <- 2     # No. cross-fitting folds for Monte Carlo analysis
KF_RD <- 5     # No. cross-fitting folds for analysis with real data
CALP  <- 0.30  # Share of sample used only for calibration
TRIM  <- 0.01  # Scores clipping threshold
NBOOT <- 1000   # Bootstrap replication (SE, real data)


# 2 Data----

data(pension)
cov_names <- c("age", "inc", "educ", "fsize", "marr", "twoearn",
               "db", "pira", "hown")
X <- as.matrix(pension[, cov_names])
D <- pension$e401
Y <- pension$net_tfa / 1000
n <- nrow(X)
cat(sprintf("n = %d, part traitée = %.3f\n", n, mean(D)))

# 3 Functions----

logit <- function(p) log(p / (1 - p))
expit <- function(z) 1 / (1 + exp(-z))
clip  <- function(p, eps = TRIM) pmin(pmax(p, eps), 1 - eps)

## 3.1 Calibration ----

platt_cal <- function(p_cal, d_cal, p_new) {
  z <- logit(clip(p_cal, 1e-6))
  fit <- suppressWarnings(glm(d_cal ~ z, family = binomial))

  expit(coef(fit)[1] + coef(fit)[2] * logit(clip(p_new, 1e-6)))
}

beta_cal <- function(p_cal, d_cal, p_new) {
  z1 <- log(clip(p_cal, 1e-6)); z2 <- -log(1 - clip(p_cal, 1e-6))
  fit <- suppressWarnings(glm(d_cal ~ z1 + z2, family = binomial))

  expit(coef(fit)[1] + coef(fit)[2] * log(clip(p_new, 1e-6))
        + coef(fit)[3] * (-log(1 - clip(p_new, 1e-6))))
}

iso_cal <- function(p_cal, d_cal, p_new) {
  o <- order(p_cal); fit <- isoreg(p_cal[o], d_cal[o])
  sf <- approxfun(fit$x, fit$yf, method = "constant", rule = 2, f = 0)

  pmin(pmax(sf(p_new), 1e-4), 1 - 1e-4)
}

## 3.2 Scoring Models----

fit_pscore <- function(kind, X, d) {
  if (kind == "logit") {
    df <- data.frame(d = d, X)
    fit <- suppressWarnings(glm(d ~ ., data = df, family = binomial))
    function(Z) clip(suppressWarnings(
      predict(fit, newdata = data.frame(Z), type = "response")), 1e-6)

  } else if (kind == "rf") {
    fit <- ranger(
      y = factor(d), x = data.frame(X), probability = TRUE,
      num.trees = 300, min.node.size = 2, num.threads = 1
    )
    function(Z) predict(fit, data = data.frame(Z))$predictions[, "1"]

  } else if (kind == "boost") {
    fit <- xgboost(
      data = X, label = factor(d), nrounds = 300, max_depth = 4,
      eta = 0.1, objective = "binary:logistic", verbose = 0,
      nthread = 1
    )
    function(Z) predict(fit, Z)

  } else stop(kind)
}

crossfit_scores <- function(X, D, K) {
  folds <- sample(rep(1:K, length.out = length(D)))
  models <- c("logit", "rf", "boost")
  out <- lapply(models, function(m)
    setNames(replicate(4, rep(NA_real_, length(D)), simplify = FALSE),
             c("raw", "platt", "beta", "iso")))
  names(out) <- models
  for (k in 1:K) {
    tr <- which(folds != k); te <- which(folds == k)
    cal_idx <- sample(tr, floor(CALP * length(tr)))
    fit_idx <- setdiff(tr, cal_idx)
    for (m in models) {
      pred  <- fit_pscore(m, X[fit_idx, , drop = FALSE], D[fit_idx])
      p_cal <- pred(X[cal_idx, , drop = FALSE])
      p_te  <- pred(X[te, , drop = FALSE])
      out[[m]]$raw[te]   <- p_te
      out[[m]]$platt[te] <- platt_cal(p_cal, D[cal_idx], p_te)
      out[[m]]$beta[te]  <- beta_cal(p_cal, D[cal_idx], p_te)
      out[[m]]$iso[te]   <- iso_cal(p_cal, D[cal_idx], p_te)

    }
  }

  out
}

crossfit_outcome <- function(X, Y, D, K) {
  folds <- sample(rep(1:K, length.out = length(D)))
  g0 <- g1 <- rep(NA_real_, length(D))
  for (k in 1:K) {
    tr <- which(folds != k); te <- which(folds == k)
    for (d in 0:1) {
      idx <- tr[D[tr] == d]
      fit <- ranger(y = Y[idx], x = data.frame(X[idx, , drop = FALSE]),
                    num.trees = 300, min.node.size = 5)
      pr <- predict(fit, data = data.frame(X[te, , drop = FALSE]))$predictions
      if (d == 0) g0[te] <- pr else g1[te] <- pr
    }

  }

  list(g0 = g0, g1 = g1)
}

## 3.3 Estimators----

estimators <- function(Y, D, p, g0, g1) {
  p <- clip(p)
  ipw  <- mean(Y * (D / p - (1 - D) / (1 - p)))
  w1   <- (D / p) / sum(D / p)
  w0   <- ((1 - D) / (1 - p)) / sum((1 - D) / (1 - p))
  ipwn <- sum(Y * w1) - sum(Y * w0)
  aipw <- mean(g1 - g0 + D * (Y - g1) / p - (1 - D) * (Y - g0) / (1 - p))
  aipwn <- mean(g1 - g0) + sum(D * (Y - g1) / p) / sum(D / p) -
    sum((1 - D) * (Y - g0) / (1 - p)) / sum((1 - D) / (1 - p))

  c(IPW = ipw, `IPW-N` = ipwn, AIPW = aipw, `AIPW-N` = aipwn)
}

## 3.4 Quality of Fit----

#' Expected Calibration Error
#' @description
#' Computes the expected calibration error (ECE) between observed binary
#' outcomes and predicted scores (probabilities) using a binned approximation:
#' \deqn{\mathrm{ECE} = \sum_{b=1}^{B} w_b \, \lvert \mathrm{acc}_b -
#' \mathrm{conf}_b \rvert}
#' where \eqn{\mathrm{acc}_b} is the mean of `d` in bin \eqn{b},
#' \eqn{\mathrm{conf}_b} is the mean of `scores` in bin \eqn{b}, and
#' \eqn{w_b} is the proportion of observations in bin \eqn{b}.
#'
#' @param d Vector of binary treatments.
#' @param scores Vector of scores.
#' @param n_bins Number of bins (default to 10).
#' @param binning Character string specifying the binning rule (either `"width"
#'  or `"quantile"`, see details).
#'
#' @details
#' The binning strategies:
#' * `"width`: equal-width bins over the interval \eqn{[0, 1]}.
#' * `"quantile`": quantile-based breaks computed from `scores`. Repeated
#'     quantiles can yield duplicate breaks; duplicates are removed
#'
#' @returns A single non-negative numeric value giving the expected calibration
#'  error. Values closer to 0 indicate better calibration.
compute_ece <- function(d,
                        scores,
                        n_bins = 10,
                        binning = c("width", "quantile")) {

  binning <- match.arg(binning)
  n <- length(d)

  # Define bin breaks
  if (binning == "width") {
    breaks <- seq(0, 1, length.out = n_bins + 1)
  } else {
    breaks <- quantile(
      scores,
      probs = seq(0, 1, length.out = n_bins + 1),
      names = FALSE
      #type = 7
    )
    breaks <- unique(breaks)
  }

  # Assign bins
  bins <- cut(
    scores,
    breaks = breaks,
    include.lowest = TRUE,
    right = TRUE
  )

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

# 4 Simulations----

## 4.1 Semi-synthetic Experiment: known ATE----

# n = 2000, p = 200, known ATE

Z <- scale(X)
z_age <- Z[, "age"]; z_inc <- Z[, "inc"]; z_edu <- Z[, "educ"]
z_fs  <- Z[, "fsize"]
marr <- X[, "marr"]; pira <- X[, "pira"]; hown <- X[, "hown"]

lin    <- -0.3 + 0.8 * z_inc + 0.5 * z_age - 0.4 * z_age^2 + 0.6 * pira +
  0.5 * z_inc * marr - 0.4 * z_fs + 0.3 * hown * z_edu
e_true <- pmin(pmax(expit(lin), 0.02), 0.98)
mu0    <- 10 + 3 * z_inc + 2 * z_age^2 + 2 * hown + 1.5 * z_edu * z_inc
tau    <- 5 + 2 * z_inc + 1.5 * pira
ATE    <- mean(tau)
cat(sprintf("ATE vrai = %.3f\n", ATE))

# Noisy proxies
set.seed(77)
jj <- sample(ncol(X), P_ADD, replace = TRUE)
X_aug <- cbind(X, Z[, jj] + matrix(rnorm(n * P_ADD), n, P_ADD))


library(future)
library(future.apply)
library(progressr)
handlers(global = TRUE)
handlers("cli")
n_cores <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_cores)

res_list <- with_progress({

  # Create a progress step counter
  p <- progressor(steps = NREP)

  future_lapply(1:NREP, function(r) {
    set.seed(1000 + r)
    sub <- sample(n, NSUB)
    Xa <- X_aug[sub, ]; e_s <- e_true[sub]
    Dr <- as.integer(runif(NSUB) < e_s)
    Yr <- mu0[sub] + tau[sub] * Dr + rnorm(NSUB, 0, 5)

    # ggplot(
    #   data = tibble(e_s = e_s, Dr = as.character(Dr)),
    #   mapping = aes(x = e_s)
    # ) +
    #   geom_histogram(mapping = aes(fill = Dr), position = "identity", alpha = .6)

    sc <- crossfit_scores(Xa, Dr, KF_MC)
    gg <- crossfit_outcome(Xa, Yr, Dr, KF_MC)
    r2g0 <- 1 - mean((mu0[sub] - gg$g0)^2) / var(mu0[sub])

    iter_rows <- list()
    for (m in names(sc)) {
      for (cb in names(sc[[m]])) {
        est <- estimators(Yr, Dr, sc[[m]][[cb]], gg$g0, gg$g1)
        iter_rows[[length(iter_rows) + 1]] <- data.frame(
          rep = r, model = m, calib = cb, estimator = names(est),
          estimate = as.numeric(est),
          ece = compute_ece(scores = sc[[m]][[cb]], d = Dr, binning = "quantile"),
          r2g0 = r2g0, ATE_true = ATE
        )
      }
    }

    p() # progress after each iteration finishes

    do.call(rbind, iter_rows)
  }, future.seed = TRUE)
})

# Combine and save results
results_mc <- list_rbind(res_list)

plan(sequential) # Close background workers

# Export results
write_csv(results_mc, paste0(output_dir, "results_mc_pension.csv"))


## 4.2 Real Data----

# Complete dataset e401

set.seed(20260712)
sc <- crossfit_scores(X, D, KF_RD)
gg <- crossfit_outcome(X, Y, D, KF_RD)

aipw_if_se <- function(Y, D, p, g0, g1) {
  p <- clip(p)
  psi <- g1 - g0 + D * (Y - g1) / p - (1 - D) * (Y - g0) / (1 - p)
  c(mean(psi), sd(psi) / sqrt(length(Y)))
}

boot_se <- function(stat_fun, B = NBOOT) {
  sd(replicate(B, stat_fun(sample.int(n, n, replace = TRUE))))
}

rows <- list()
for (m in names(sc)) for (cb in c("raw", "beta")) {
  p <- clip(sc[[m]][[cb]])
  est <- estimators(Y, D, p, gg$g0, gg$g1)
  ifr <- aipw_if_se(Y, D, p, gg$g0, gg$g1)

  # Calibration error (quantile-binned ECE) for this model / calibration branch
  ece <- compute_ece(scores = sc[[m]][[cb]], d = Dr, binning = "quantile")
  ece_se <- boot_se(function(i)
    compute_ece(scores = sc[[m]][[cb]][i], d = Dr[i], binning = "quantile"))


  se <- c(
    IPW = boot_se(function(i)
      mean(Y[i] * (D[i] / p[i] - (1 - D[i]) / (1 - p[i])))),
    `IPW-N` = boot_se(function(i) {
      w1 <- (D[i] / p[i]) / sum(D[i] / p[i])
      w0 <- ((1 - D[i]) / (1 - p[i])) / sum((1 - D[i]) / (1 - p[i]))
      sum(Y[i] * w1) - sum(Y[i] * w0) }),
    AIPW = ifr[2],
    `AIPW-N` = boot_se(function(i)
      mean(gg$g1[i] - gg$g0[i]) +
        sum(D[i] * (Y[i] - gg$g1[i]) / p[i]) / sum(D[i] / p[i]) -
        sum((1 - D[i]) * (Y[i] - gg$g0[i]) / (1 - p[i])) /
        sum((1 - D[i]) / (1 - p[i])))
  )

  for (k in names(est)) {
    rows[[length(rows) + 1]] <- data.frame(
      model = m, calib = cb,
      estimator = ifelse(k == "AIPW-N", "AIPW-H", k),
      est = est[k], se = se[k],
      ece = ece, ece_se = ece_se
    )
  }
}
real <- do.call(rbind, rows) |> tibble::as_tibble()

# Export results
write_csv(real, file = paste0(output_dir, "estimates_real_pension.csv"))


# 5 Results----

## 5.1 Semi-synthetic Experiment----

results_mc <- read_csv(paste0(output_dir, "results_mc_pension.csv"))

### Figure 2----

true_ate <- results_mc$ATE_true[1]

tb_plot <- results_mc |>
  filter(
    str_detect(estimator, "^AIPW"),
    calib %in% c("raw", "beta"),
    (calib != "beta" | estimator != "AIPW-N")
  ) |>
  mutate(
    lab = case_when(
      calib == "raw" & estimator == "AIPW" ~ "Uncalibrated",
      calib == "raw" & estimator == "AIPW-N" ~ "Hájek-normalized",
      calib == "beta" & estimator == "AIPW" ~ "Beta-calibrated",
      TRUE ~ "Error"
    ),
    lab = factor(lab, levels = c("Uncalibrated", "Hájek-normalized", "Beta-calibrated", "Error")),
    model = factor(
      model,
      levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random Forest", "Boosting")
    )
  )

p_bp_semisynth_mc <-
  ggplot(
  data = tb_plot,
  mapping = aes(x = model, y = estimate, fill = lab)
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.25,
    position = position_dodge(width = 0.75)
  ) +
  geom_boxplot(
    outliers = FALSE,
    width = 0.65,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(
    mapping = aes(yintercept = true_ate, linetype = "True ATE"),
    colour = "#009E73", linewidth = 1.5
  ) +
  labs(
    x = NULL, y = "Estimated ATE (AIPW)"
  ) +
  scale_fill_manual(
    NULL,
    values = c(
      "Uncalibrated" = "darkgrey",
      "Hájek-normalized" = "#0072B2",
      "Beta-calibrated" = "#D55E00",
      "Error" = "gray")
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("True ATE" = "dashed")
  ) +
  theme(text = element_text(size = 16)) +
  theme_paper() +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major.y = element_line(colour = "gray", linewidth = .2)
  )

if (isTRUE(export_graphs)) {
  ggplot2_to_pdf(
    plot = p_bp_semisynth_mc, path = fig_dir, filename = "semisynth_mc",
    width = 8, height = 3
  )
}

### Table 4. descriptive statistics on the MC simulations----

tb_sim_semisync_res <- results_mc |>
  group_by(model, calib, estimator) |>
  mutate(bias = estimate - ATE_true) |>
  summarise(
    ate_mean = mean(estimate),
    ate_sd = sd(estimate),
    ece_mean = mean(ece),
    ece_sd = sd(ece),
    r2g0_mean = mean(r2g0),
    r2g0_sd = sd(r2g0),
    bias_mean = mean(bias),
    bias_sd = sd(bias)
  ) |>
  mutate(
    model = factor(
      model,
      levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random Forest", "Boosting")
    ),
    calib = factor(
      calib,
      levels = c("raw", "platt", "beta", "iso"),
      labels = c("Uncalibrated", "Platt", "Beta", "Isotonic")
    ),
    estimator = factor(
      estimator,
      levels = c("IPW", "IPW-N", "AIPW", "AIPW-N"),
      labels = c("IPW", "IPW (Hájek-normalized)", "AIPW", "AIPW (Hájek-normalized)")
    )
  ) |>
  arrange(model, calib, estimator)

#' Format numbers
#'
#' |x| >= 100: 0 digit;
#' |x| >= 10: 1 digit;
#' otherwise: 3 digits
fmt_num <- function(x) {
  dplyr::case_when(
    abs(x) >= 100 ~ sprintf("%.0f", x),
    abs(x) >= 10  ~ sprintf("%.1f", x),
    TRUE          ~ sprintf("%.3f", x)
  )
}

tb_to_save <- tb_sim_semisync_res |>
  ungroup() |>
  dplyr::select(-r2g0_mean, -r2g0_sd) |>
  filter(
    estimator %in% c("AIPW", "AIPW (Hájek-normalized)"),
    calib %in% c("Uncalibrated", "Beta")
    ) |>
  filter(
    !(calib == "Beta" & estimator == "AIPW (Hájek-normalized)")
  ) |>
  mutate(
    across(where(is.numeric), fmt_num),
    across(matches("_sd$"), ~ paste0("(", .x, ")"))
  ) |>
  mutate(
    ate = str_c(ate_mean, ate_sd, sep = " "),
    ece = str_c(ece_mean, ece_sd, sep = " ")
  ) |>
  dplyr::select(-ate_mean, -ate_sd, -ece_mean, -ece_sd) |>
  dplyr::select(-bias_mean, -bias_sd) |>
  mutate(
    lab = case_when(
      calib == "Uncalibrated" & estimator == "AIPW" ~ "Uncalibrated",
      calib == "Uncalibrated" & estimator == "AIPW (Hájek-normalized)" ~ "Hájek-normalized",
      calib == "Beta" & estimator == "AIPW" ~ "Beta-calibrated",
      TRUE ~ "Error"
    ),
    lab = factor(lab, levels = c("Uncalibrated", "Hájek-normalized", "Beta-calibrated", "Error"))
  ) |>
  relocate(lab, .after = "model") |>
  dplyr::select(-calib, -estimator)


table_4_lines <- apply(tb_to_save, 1, function(x) {
  paste(paste(x, collapse = " & "), "\\\\")
})

writeLines(table_4_lines, file.path(tbl_dir, "table4-body.tex"))

  # kableExtra::kbl(digits = 2) |>
  # kableExtra::kable_styling(
  #   bootstrap_options = c("striped", "hover", "condensed"),
  #   latex_options = c("striped", "hold_position")
  # )

### Figure F.1 (Appendix)----

y_lim <- 20

tb_plot <- results_mc |>
  filter(
    str_detect(estimator, "^IPW"),
    calib %in% c("raw", "beta")
  ) |>
  mutate(
    calib_old = calib,
    model = factor(
      model,
      levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random Forest", "Boosting")
    ),
    estimator = factor(
      estimator,
      levels = c("IPW", "IPW-N"),
      labels = c("IPW", "Normalized IPW (Hájek)")
    )
  ) |>
  group_by(estimator, model, calib) |>
  mutate(
    estimate_mean = mean(estimate, na.rm = TRUE),
    # Replace estimates > y_lim with NA so ggplot keeps the factor level structure
    calib = ifelse(estimate_mean > y_lim, "to-remove", calib)
  ) |>
  ungroup() |>
  mutate(
    calib = factor(
      calib,
      levels = c("to-remove", "raw", "beta"),
      labels = c("to-remove", "Uncalibrated", "Beta-calibrated")
    )
  )


tb_arrows <- tb_plot |>
  filter(estimate_mean > !!y_lim) |>
  distinct(estimator, model, calib_old, calib, estimate_mean) |>
  mutate(
    calib = factor(
      calib_old,
      levels = c("raw", "beta"),
      labels = c("Uncalibrated", "Beta-calibrated")
    )
  ) |>
  mutate(
    x = case_when(
      calib == "Uncalibrated"    ~ as.numeric(model) - 0.125,
      calib == "Beta-calibrated" ~ as.numeric(model) + 0.125
    )
  )
tb_arrows

colour_arrow <- "#505050"

p_bp_semisynth_ipw <-
  ggplot(
  data = tb_plot,
  mapping = aes(x = model, y = estimate, fill = calib, colour = calib)
) +
  stat_boxplot(
    geom = "errorbar",
    width = 0.5
  ) +
  geom_boxplot(
    outliers = FALSE,
    width = 0.5
  ) +
  # # Upward arrows pointing to the top edge (y = 17.5 to 19.8)
  geom_segment(
    data = tb_arrows,
    mapping = aes(
      x = x,
      xend = x,
      y = y_lim - 2.5,
      yend = y_lim - 0.2,
      group = calib
    ),
    arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
    linewidth = 0.8,
    inherit.aes = FALSE,
    show.legend = FALSE,
    colour = colour_arrow
  ) +
  # Text labels displaying the average ATE below the arrow
  geom_text(
    data = tb_arrows,
    mapping = aes(
      x = x, y = y_lim - 3.5,
      label = sprintf("%.1f", estimate_mean),
      group = calib
    ),
    inherit.aes = FALSE,
    position = position_dodge(width = 0.75),
    size = 3, fontface = "bold",
    show.legend = FALSE,
    colour = colour_arrow
  ) +
  geom_hline(
    mapping = aes(yintercept = true_ate, linetype = "True ATE"),
    colour = "#009E73", linewidth = 1
  ) +
  facet_wrap(~estimator) +
  coord_cartesian(ylim = c(0, y_lim)) + # Use coord_cartesian instead of ylim() to avoid clipping warnings
  labs(
    x = NULL, y = "Estimated ATE"
  ) +
  scale_fill_manual(
    NULL,
    values = c(
      "Uncalibrated" = "darkgrey",
      "Beta-calibrated" = "#D55E00"
    ),
    drop = FALSE # Keep all legend categories even if empty
  ) +
  scale_colour_manual(
    values = c(
      "Uncalibrated" = "#0A0A0A",
      "Beta-calibrated" = "#0A0A0A",
      "to-remove" = NA
    ),
    guide = "none"
  ) +
  scale_linetype_manual(
    name = NULL,
    values = c("True ATE" = "dashed")
  ) +
  theme_paper() +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major.y = element_line(colour = "gray", linewidth = .2)
  )

if (isTRUE(export_graphs)) {
  ggplot2_to_pdf(
    plot = p_bp_semisynth_ipw, path = fig_dir, filename = "semisynth_ipw",
    width = 9, height = 4
  )
}


### Supplementary figure----
p_all <-
  ggplot(
    data = results_mc |>
      mutate(
        calib = factor(calib, levels = c("raw", "platt", "beta", "iso")),
        model = factor(model, levels = c("logit", "rf", "boost"),
                       labels = c("Logistic", "Random Forest", "Boosting")),
        estimator = factor(estimator, levels = c("IPW", "IPW-N", "AIPW", "AIPW-N"))
      ),
    mapping = aes(x = calib, y = estimate, fill = estimator)
  ) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  geom_hline(yintercept = true_ate, linetype = "dashed", colour = "black", linewidth = 0.7) +
  facet_grid(estimator ~ model, scales = "free_y") +
  labs(
    title = "ATE Estimates Across Models, Estimators, and Calibration Schemes",
    x = "Calibration Scheme",
    y = "Estimate",
    fill = "Estimator"
  ) +
  theme_paper()


## 5.2 Real Data----

real <- read_csv(paste0(output_dir, "estimates_real_pension.csv"))

### Table 5 ----
tb_401k_aipw <- real |> filter(str_detect(estimator, "^AIPW")) |>
  mutate(
    ate = str_c(fmt_num(est), " (", fmt_num(se), ")"),
    ece = str_c(fmt_num(ece), " (", fmt_num(ece_se), ")"),
  ) |>
  dplyr::select(-est, -se, -ece_se) |>
  mutate(
    type = case_when(
      calib == "raw" & estimator == "AIPW" ~ "Uncalibrated",
      calib == "raw" & estimator == "AIPW-H" ~ "Hájek",
      calib == "beta" & estimator == "AIPW" ~ "Beta-calibrated",
      TRUE ~ NA_character_
    ),
    model = factor(
      model, levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random forest", "Boosting")
    )
  ) |>
  filter(!is.na(type)) |>
  dplyr::select(-estimator, -calib) |>
  mutate(outcome = "Random Forest") |>
  dplyr::select(outcome, model, type, ate, ece)


table_5_lines <- apply(tb_401k_aipw, 1, function(x) {
  paste(paste(x, collapse = " & "), "\\\\")
})

writeLines(table_5_lines, file.path(tbl_dir, "table5-body.tex"))


### Figure 3----

tb_plot <- real |>
  filter(
    str_detect(estimator, "^AIPW"),
    calib %in% c("raw", "beta"),
    (calib != "beta" | estimator != "AIPW-H")
  ) |>
  mutate(
    lab = case_when(
      calib == "raw" & estimator == "AIPW" ~ "Uncalibrated",
      calib == "raw" & estimator == "AIPW-H" ~ "Hájek-normalized",
      calib == "beta" & estimator == "AIPW" ~ "Beta-calibrated",
      TRUE ~ "Error"
    ),
    lab = factor(lab, levels = c("Uncalibrated", "Hájek-normalized", "Beta-calibrated", "Error")),
    model = factor(
      model,
      levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random Forest", "Boosting")
    )
  )

p_401k_real <- ggplot(
  data = tb_plot,
  mapping = aes(x = est, y = fct_rev(model), colour = lab)
) +
  geom_errorbar(
    mapping = aes(xmin = est - se, xmax = est + se),
    width = .15,
    position = position_dodge(width = .4)
  ) +
  geom_point(position = position_dodge(width = .4)) +
  labs(
    x = NULL, y = "Estimated ATE (AIPW)"
  ) +
  scale_colour_manual(
    NULL,
    values = c(
      "Uncalibrated" = "darkgrey",
      "Hájek-normalized" = "#0072B2",
      "Beta-calibrated" = "#D55E00",
      "Error" = "gray")
  ) +
  theme_paper() +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major.x = element_line(colour = "gray", linewidth = .2)
  )


if (isTRUE(export_graphs)) {
  ggplot2_to_pdf(
    plot = p_401k_real, path = fig_dir, filename = "401k_real",
    width = 7, height = 3
  )
}

### Table with descriptive statistics on real data----

real |>
  mutate(
    model = factor(
      model,
      levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random Forest", "Boosting")
    ),
    calib = factor(
      calib,
      levels = c("raw", "platt", "beta", "iso"),
      labels = c("Uncalibrated", "Platt", "Beta", "Isotonic")
    ),
    estimator = factor(
      estimator,
      levels = c("IPW", "IPW-N", "AIPW", "AIPW-H"),
      labels = c("IPW", "IPW (Hájek-normalized)", "AIPW", "AIPW (Hájek-normalized)")
    )
  ) |>
  arrange(model, calib, estimator)


  # kableExtra::kbl(digits = 2) |>
  # kableExtra::kable_styling(
  #   bootstrap_options = c("striped", "hover", "condensed"),
  #   latex_options = c("striped", "hold_position")
  # )

### Figure F.2 (Appendix)----

tb_plot <- real |>
  filter(
    str_detect(estimator, "^IPW"),
    calib %in% c("raw", "beta")
  ) |>
  mutate(
    model = factor(
      model,
      levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random Forest", "Boosting")
    ),
    estimator = factor(
      estimator,
      levels = c("IPW", "IPW-N"),
      labels = c("IPW", "Normalized IPW (Hájek)")
    ),
    calib = factor(
      calib,
      levels = c("raw", "beta"),
      labels = c("Uncalibrated", "Beta-calibrated")
    )
  )

p_401k_ipw <- ggplot(
  data = tb_plot,
  mapping = aes(x = est, y = fct_rev(model), colour = calib)
) +
  geom_errorbar(
    mapping = aes(xmin = est - se, xmax = est + se),
    width = .15,
    position = position_dodge(width = .4)
  ) +
  geom_point(position = position_dodge(width = .4)) +
  geom_vline(xintercept = 0) +
  facet_wrap(~estimator) +
  labs(
    x = "ATE (thousands of dollars)", y = NULL
  ) +
  scale_colour_manual(
    NULL,
    values = c(
      "Uncalibrated" = "darkgrey",
      "Beta-calibrated" = "#D55E00"
    ),
    drop = FALSE # Keep all legend categories even if empty
  ) +
  theme_paper() +
  theme(
    panel.background = element_rect(fill = NA),
    panel.grid.major.x = element_line(colour = "gray", linewidth = .2)
  )

if (isTRUE(export_graphs)) {
  ggplot2_to_pdf(
    plot = p_401k_ipw, path = fig_dir, filename = "401k_ipw",
    width = 9, height = 4
  )
}

### Figure E.1 (Appendix)----

reliability_diagram <- function(p, d, num_bins = 10) {
  # Generate quantile probability breaks (0, 0.1, 0.2, ..., 1.0 for 10 bins)
  breaks_q <- quantile(p, probs = seq(0, 1, length.out = num_bins + 1), na.rm = TRUE)

  # Ensure strictly unique breaks if duplicate predictions cause identical quantiles
  breaks_q <- unique(breaks_q)

  data.frame(y_pred = p, y_obs = d) |>
    mutate(
      # Bin by empirical quantiles instead of equal width
      bin = cut(y_pred, breaks = breaks_q, include.lowest = TRUE)
    ) |>
    group_by(bin) |>
    summarise(
      n_obs     = n(),
      mean_pred = mean(y_pred),
      mean_obs  = mean(y_obs),
      # Standard error for proportion
      se_obs    = sqrt((mean_obs * (1 - mean_obs)) / n_obs),
      .groups   = "drop"
    )
}

calib_curve <- list()
for (m in names(sc)) {
  for (cal in names(sc[[m]])) {
    tmp_calib_curve <- reliability_diagram(sc[[m]][[cal]], D)
    calib_curve[[length(calib_curve) + 1]] <- tmp_calib_curve |>
      mutate(model = m, calib = cal)
  }
}

calib_curve <- calib_curve |> list_rbind()

tb_plot <- calib_curve |>
  filter(calib %in% c("raw", "beta")) |>
  mutate(
    model = factor(
      model,
      levels = c("logit", "rf", "boost"),
      labels = c("Logit", "Random Forest", "Boosting")
    ),
    calib = factor(
      calib,
      levels = c("raw", "beta"),
      labels = c("Uncalibrated", "Beta-calibrated")
    )
  )

p_401k_calibration <- ggplot(
  data = tb_plot,
  mapping = aes(x = mean_pred, y = mean_obs, colour = calib)
) +
  geom_line() +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~model) +
  labs(
    x = "Predicted scores", y = "Observed frequency of D=1"
  ) +
  scale_colour_manual(
    NULL,
    values = c(
      "Uncalibrated" = "darkgrey",
      "Beta-calibrated" = "#D55E00"
    ),
    drop = FALSE # Keep all legend categories even if empty
  ) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  theme_paper()


if (isTRUE(export_graphs)) {
  ggplot2_to_pdf(
    plot = p_401k_calibration, path = fig_dir, filename = "401k_calibration",
    width = 9, height = 4
  )
}
