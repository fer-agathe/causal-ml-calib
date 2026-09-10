# Propensity Score Calibration in Causal Machine Learning: What It Fixes, and What It Does Not

[![R](https://img.shields.io/badge/Language-R-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Quarto](https://img.shields.io/badge/Built%20with-Quarto-1a9c8a?logo=quarto&logoColor=white)](https://quarto.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

Replication code and companion e-book for:

> Fernandes Machado, A., Charpentier, A., Flachaire, E., and Gallic, E. *Propensity Score Calibration in Causal Machine Learning: What It Fixes, and What It Does Not.*

**Abstract.** A recent literature shows that recalibrating machine learning estimates of the propensity score improves double/debiased machine learning (DML) estimators of average treatment effects. This paper asks what calibration can and cannot fix, and how it compares with the weight-stabilization devices routinely used by applied econometricians. Three answers emerge from an elementary bias decomposition and a factorial Monte Carlo design. First, calibration matters for the bias only when both nuisance functions are poorly estimated (with an accurate outcome regression, double robustness protects the AIPW estimator against arbitrarily miscalibrated scores), yet it reduces the standard deviation when the propensity learner produces extreme scores. Second, calibration dominates common practice under limited overlap: Hájek normalization tames the variance explosions caused by extreme scores but leaves most of the bias, whereas calibration removes both at once and without a tuning parameter, unlike trimming. Third, calibration is necessary but not sufficient for uniform unbiasedness over outcome functions: the bias of weighting estimators built on a calibrated score equals an average within-stratum covariance between score errors and (residualized) potential outcomes, so a perfectly calibrated score can still deliver a first-order biased estimator.

**Keywords:** Average treatment effect; Double machine learning; Propensity score; Calibration; AIPW; Overlap.

**JEL codes:** C14, C21, C52.

A companion e-book with detailed walkthroughs of the simulation design, code, and results is available at: **[https://fer-agathe.github.io/causal-ml-calib/](https://fer-agathe.github.io/causal-ml-calib/)** (rendered from `ebook/`, deployed to `docs/`).

---

## Repository structure

```
.
├── ebook/            # Companion e-book (Quarto) walking through the simulation design and results
├── docs/             # Rendered e-book (GitHub Pages)
│
├── scripts/
│   ├── simulation/   # Modular Monte Carlo pipeline (DGP → learners → calibration →
│   │                 # metrics → estimators), run via simulation.R, plus figs-and-tbls.R
│   │                 # to build the paper's figures/tables
│   └── 401k/         # Empirical application (401(k) eligibility)
│
├── output/           # Simulation results, generated figures (PDF) and tables (.tex)
│
├── main_v3_2_.tex    # Paper source (LaTeX)
└── LICENSE
```

## Reproducing the results

### Requirements

- **R** (≥ 4.2 recommended), with packages including `tidyverse`, `MASS`, `sandwich`, `philentropy`, `gmish`, `ebal`, `pbapply`, `tikzDevice`.
- **Quarto** (to render the companion e-book in `ebook/`).
- A LaTeX distribution (to render `tikzDevice` figures).

### 1. Monte Carlo simulation

From the repository root:

```r
source("scripts/simulation/simulation.R")
```

This sources the modular pipeline (`01_dgp.R` → `07_simulation_helper.R`) and runs the factorial Monte Carlo design, writing raw results to `output/simul-klassen/`.

Then generate the paper's figures and tables:

```r
source("scripts/simulation/figs-and-tbls.R")
```

This produces the figures in `output/figs/` (PDF) and the table bodies in `output/tbl/`.

### 2. Empirical application (401(k) eligibility)

```r
source("scripts/401k/application_pension.R")
```

Writes estimates to `output/401k/`.


## Paper

The working paper will be abailable soon.

## Authors

- Agathe Fernandes Machado — Université du Québec à Montréal
- Arthur Charpentier — Université du Québec à Montréal
- Emmanuel Flachaire — Aix-Marseille Univ., CNRS, AMSE
- Ewen Gallic — Aix-Marseille Univ., CNRS, AMSE

## Funding
Agathe Fernandes Machado acknowledges funding from OBVIA, Aix-Marseille School of Economics, Université du Québec à Montréal and the Centre de Recherches Mathématiques. Emmanuel Flachaire and Ewen Gallic acknowledge funding from the French government under the "France 2030" investment plan managed by the French National Research Agency (reference: ANR-17-EURE-0020), and from the Excellence Initiative of Aix-Marseille University – A*MIDEX.

## License

This repository is licensed under the **GNU General Public License v3.0 (GPL-3.0)** — see [LICENSE](LICENSE) for the full text. In short: you are free to use, study, modify, and redistribute this code, including for derivative work, as long as derivative work is also released under GPL-3.0 with source code available.
