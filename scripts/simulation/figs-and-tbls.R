# Figures and tables for the Monte Carlo simulations

library(tidyverse)
library(tikzDevice)
source("scripts/theme.R")

output_dir <- "output/simul-klassen/" # Simulation results
output_dir_figs <- "output/figs/"     # Figures
output_dir_tbl <- "output/tbl"        # Tables

if (!dir.exists(output_dir_figs)) dir.create(output_dir_figs, recursive = TRUE)
if (!dir.exists(output_dir_tbl)) dir.create(output_dir_tbl, recursive = TRUE)

# Load Results----

cat("Processing and assembling results structure...\n")

res <- list()

for (i_grid in 1:6) {
  cat(
    sprintf(
      "Load results from simulation (%d/%d)...\n\n",
      i_grid, 6
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


## Functions for tables----

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


# Combine into master dataframes
res_ate     <- map(res, "ate") |> list_rbind()
res_gate    <- map(res, "gate") |> list_rbind()
res_metrics <- map(res, "metrics") |> list_rbind()

# Figures and Tables----

## Figure 1----


p_hist <- ggplot(
  data = res_ate |>
    filter(
      g_name == "lgbm", m_name == "lgbm",
      calib_name %in% c("uncalibrated", "isotonic"),
      n == 2000, p == 200, Rd2 == 0.2
    ) |>
    mutate(
      delta_ate_aipw = ate_aipw - ate_true,
      calib_name = factor(
        calib_name,
        levels = c("uncalibrated", "isotonic"),
        labels = c("Without calibration", "With calibration (isotonic)")
      )
    ),
  mapping = aes(x = delta_ate_aipw)
) +
  geom_histogram(
    colour = "black", fill = "gray", linewidth = .4
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ calib_name) +
  coord_cartesian(xlim = c(-1, .25)) +
  scale_x_continuous(breaks = seq(-1, .2, by = .2)) +
  labs(
    x = "$\\hat{\\tau}_{\\text{aipw}} - \\hat{\\tau}$",
    y = "Frequency"
  ) +
  theme_paper() +
  theme(
    base_size = 14,
    panel.grid.major.y = element_line(colour = "gray", linewidth = .2)
    )

if (!dir.exists(output_dir_figs)) dir.create(output_dir_figs, recursive = TRUE)
# ggsave(p_hist, file = "output/figs/fig1-histogram.pdf", width = 9, height = 4)

ggplot2_to_pdf(
  plot = p_hist, path = output_dir_figs, filename = "fig1-histogram",
  height = 3.5, width = 9, keep_tex = FALSE
)

## Table 1----
# p = 200, R2d = 0.2

g_levels <- c(ols = "OLS", lasso = "Lasso", rf = "RF", lgbm = "LGBM")
m_levels <- c(
  logit = "Logit", rf_classif = "RF-Classif", rf_reg = "RF-Reg", lgbm = "LGBM"
)

#' Format one cell (Bias (SD)) for a given configuration
cell <- function(data, estimator_col, calib) {
  e <- data |>
    filter(calib_name == calib) |>
    mutate(err = .data[[estimator_col]] - ate_true) |>
    pull(err)
  sprintf("%s & (%s)", fmt_num(mean(e)), fmt_num(sd(e)))
}

#' Body of the table
#'
#' @param p_val Number of variables in the simulation
#' @param rd2_val Numeric value specifying the input \exp{R^2_d}.
#' @param spec: Named list with two elements: `col`, the column with the
#'   estimate, and `calib`, the name of the calibration technique.
make_body <- function(p_val,
                      rd2_val,
                      spec) {
  sub <- res_ate |> filter(p == p_val, Rd2 == rd2_val)
  lines <- c()
  for (gi in seq_along(g_levels)) {
    g <- names(g_levels)[gi]
    for (mi in seq_along(m_levels)) {
      m <- names(m_levels)[mi]
      s <- sub |> filter(g_name == g, m_name == m)
      cells <- vapply(spec, function(sp) cell(s, sp$col, sp$calib), character(1))
      lead <- if (mi == 1) g_levels[gi] else ""
      lines <- c(lines, paste0(lead, " & ", m_levels[mi], " & ",
                               paste(cells, collapse = " & "), " \\\\"))
    }
    if (gi < length(g_levels)) lines <- c(lines, "\\addlinespace")
  }
  lines
}

#' Write lines in a tex file
write_body <- function(lines, out_dir, file) {
  writeLines(lines, file.path(out_dir, file))
  cat("==", file, "==\n")
  cat(lines, sep = "\n")
  cat("\n")
}

spec_t1 <- list(
  uncal = list(col = "ate_aipw", calib = "uncalibrated"),
  platt = list(col = "ate_aipw", calib = "platt"),
  beta  = list(col = "ate_aipw", calib = "beta"),
  iso   = list(col = "ate_aipw", calib = "isotonic")
)


write_body(
  lines = make_body(200, 0.2, spec_t1),
  out_dir = output_dir_tbl,
  file = "table1-body.tex"
)

## Table 2----
# p = 200, R2d = 0.8 (limited overlap)

spec_t2 <- list(
  uncal   = list(col = "ate_aipw",      calib = "uncalibrated"),
  hajek   = list(col = "ate_aipw_norm", calib = "uncalibrated"),
  clipped = list(col = "ate_aipw",      calib = "uncalib-platt-clip"),
  iso     = list(col = "ate_aipw",      calib = "isotonic")
)
write_body(
  lines = make_body(200, 0.8, spec_t2),
  out_dir = output_dir_tbl,
  file = "table2-body.tex"
)

## Table 3----


table_3 <- res_ate |>
  filter(
    n == 2000, p == 200,
    calib_name %in% c("isotonic", "joint-platt")
  ) |>
  filter(
    (m_name == "lgbm" & g_name == "lgbm" & Rd2 == .8) |
      (m_name == "rf_classif" & g_name == "lgbm" & Rd2 == .8) |
      (m_name == "lgbm" & g_name == "lgbm" & Rd2 == .2)
  ) |>
  mutate(
    m_name = factor(m_name, levels = c("lgbm", "rf_classif"), labels = c("LGBM", "RF-Classif")),
    g_name = factor(g_name, levels = c("lgbm"), labels = c("LGBM")),
    calib_name = factor(calib_name, levels = c("isotonic", "joint-platt"))
  ) |>
  group_by(g_name, m_name, Rd2, calib_name) |>
  summarise(
    bias = mean(ate_aipw - ate_true) |> fmt_num(),
    sd = sd(ate_aipw - ate_true) |> fmt_num(),
    sd = str_c("(", sd, ")"),
    .groups = "drop"
  ) |>
  arrange(desc(Rd2), m_name) |>
  pivot_wider(names_from = calib_name, values_from = c("bias", "sd"), names_vary = "slowest")

out_file <- "table3-body.tex"
table_3_lines <- apply(table_3, 1, function(x) {
  # paste(c(x, "\\\\"), collapse = " & ")

  paste(paste(x, collapse = " & "), "\\\\")
})

writeLines(table_3_lines, file.path(output_dir_tbl, out_file))


## Table C.1 (Appendix)----
# p = 20, R2d = 0.2


write_body(
  lines = make_body(20, 0.2, spec_t1),
  out_dir = output_dir_tbl,
  file = "tableC1-body.tex"
)

## Table D.1 (Appendix)----
# p = 200, R2d = 0.8

g_levels <- c(ols = "OLS", lasso = "Lasso", rf = "RF", lgbm = "LGBM")
m_levels <- c(
  logit = "Logit", rf_classif = "RF-Classif", rf_reg = "RF-Reg", lgbm = "LGBM"
)

spec_t1 <- list(
  uncal = list(col = "ate_aipw", calib = "uncalibrated"),
  platt = list(col = "ate_aipw", calib = "platt"),
  beta  = list(col = "ate_aipw", calib = "beta"),
  iso   = list(col = "ate_aipw", calib = "isotonic")
)


write_body(
  lines = make_body(200, 0.8, spec_t1),
  out_dir = output_dir_tbl,
  file = "table7-body.tex"
)


## Figure B.1 (Appendix)----

set.seed(1234)
source("scripts/simulation/01_dgp.R")

grid <- expand_grid(
  Rd2 = c(.2, .5, .8),
  Ry2 = .5
)

format_data_belloni_example <- function(x) {
  X <- as_tibble(x$X)
  colnames(X) <- str_replace(colnames(X), "^V", "X")
  tibble(
    y = x$y,
    d = x$d,
    p = x$m_0,
    X
  )
}

scale <- 1
file_name <- "fig-hist-true-scores"
out_folder <- "output/figs/"
width_tikz <- 5*scale
height_tikz <- 2*scale
tikz(paste0(out_folder, file_name, ".tex"), width = width_tikz, height = height_tikz)

range_probas <- NULL
par(mfrow = c(1,3))
for (i in 1:nrow(grid)) {
  x <- gen_data_belloni(
    n = 2000,
    p = 200,
    Rd2 = grid$Rd2[i],
    Ry2 = grid$Ry2[i],
    rho = 0.5,
    seed = 123
  )
  tb <- format_data_belloni_example(x)
  hist(
    tb |> filter(d==0) |> pull("p"),
    col = alpha("red", .5),
    main = str_c("$R_d=", grid$Rd2[i], "$, $R_y=", grid$Ry2[i], "$"), xlim = c(0,1),
    xlab = NULL,
    font.main = 1)
  hist(tb |> filter(d==1) |> pull("p"), col = alpha("blue", .5), add = TRUE)
  range_probas <- bind_rows(
    range_probas,
    tibble(
      Rd2 = grid$Rd2[i],
      min_0 = tb |> filter(d == 0) |> pull("p") |> min(),
      max_0 = tb |> filter(d == 0) |> pull("p") |> max(),
      min_1 = tb |> filter(d == 1) |> pull("p") |> min(),
      max_1 = tb |> filter(d == 1) |> pull("p") |> max()
    )
  )
}
rm(x, tb, grid, i)

dev.off()
plot_to_pdf(
  filename = file_name,
  path = out_folder, keep_tex = FALSE
)



