# Post-selection inference for linear regression (code for submission)

This folder contains R scripts used to run the simulation experiments and the application example for a linear-regression post-selection inference project.

## Files

- `source code.R`
  - Core functions used across the experiments and the application example.
  - Includes helpers for truncated normal calculations on unions of intervals and multiple post-selection p-value / CI routines.
- `application.R`
  - A small, concrete example using `fpp3::us_change` to:
    1) select a model via an information criterion (AIC/BIC/AICc),
    2) compute p-values under different approaches (selective vs naive),
    3) compute confidence intervals.
- `experiment Figure 1 6.R`, `experiment Figure 2 3 5.R`, `experiment Figure 4.R`
  - Simulation scripts used to generate figures in the paper.
- Data/helper artifacts
  - `*.RData`, `*.rds` files are intermediate saved objects used by the figure scripts.

## Requirements

These scripts are written as plain R scripts (not an installed package). You will need R plus the relevant CRAN packages.

The application and core methods rely on packages including:

- `intervals`, `nleqslv`, `leaps`, `natural`, `MASS`
- `parallel`, `foreach`, `doParallel`, `bigstatsr`
- `fpp3` (for the application dataset `us_change`)

The simulation/figure scripts additionally use packages including:

- `mvtnorm`, `ggplot2`, `tidyr`, `dplyr`, `patchwork`, `latex2exp`
- `sets`, `future`, `future.apply`

Install (example):

```r
install.packages(c(
  "fpp3",
  "intervals","nleqslv","leaps","natural","MASS",
  "parallel","foreach","doParallel","bigstatsr",
  "mvtnorm","ggplot2","tidyr","dplyr","patchwork","latex2exp",
  "sets","future","future.apply"
))
```

## Quick start: run the application example

From this folder in an R session:

```r
# Optional: define these if you want the script to be self-contained.
# (In the provided application loops, new_x is not used unless contrast1 == 'new')
# sigmaKnown <- TRUE
# new_x <- rep(NA_real_, 4)

source("application.R")
```

Output is printed to the console (p-values and intervals).

## Reproducing the simulation figures

Each figure script is intended to be run as a standalone R script:

```r
source("experiment Figure 1 6.R")
source("experiment Figure 2 3 5.R")
source("experiment Figure 4.R")
```

Depending on your machine and the number of replications, these can take time.

## Notes

- `source code.R` currently contains a small sanity-check/plotting block for the truncated normal helpers. If you only want to load functions (no plots), you may want to comment out that block.
- Some scripts refer to variables such as `sigmaKnown` and `new_x`. Where relevant, the scripts include comments explaining how they are used.
