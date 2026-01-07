################################################################################
# Experiment: Figure 1 & 6 - Coverage comparison with unknown variance
#
# Purpose:
#   Replicates simulation results from Section 4 of the paper.
#   Compares empirical coverage of confidence intervals for prediction at new
#   points under two approaches:
#     - Adjusted (selective): saturated_interval()
#     - Unadjusted (naive): naive_interval()
#
# Setup:
#   - p = 8 predictors (true model has first 3 active: beta = [1,2,3,0,...,0])
#   - n = 50 observations
#   - Independent design (diagonal correlation matrix)
#   - Model selection via AIC
#   - Variance unknown (estimated from selected model)
#
# Output:
#   - Coverage plot comparing adjusted vs unadjusted intervals across 10 new points
#   - Table of model size frequencies
#   - Coverage stratified by selected model size
#
# Note: Uses parallel computation via mclapply() for 500 replications.
################################################################################

rm(list = ls())
source('source code.R')
## --- Load required packages ---
library(intervals)
library(sets)
library(nleqslv)
library(leaps)
library(natural)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(bigstatsr)
library(mvtnorm)
library(intervals)
library(nleqslv)
library(MASS)
library(leaps)
library(future)
library(future.apply)
library(parallel)
library(ggplot2)

## --- Initialize storage for results ---
se_List = list()
c1_List = list()
c2_List = list()
ynew_List = list()

## --- Main simulation function ---
# Performs one replication:
#   1. Generate data from true model
#   2. Select model via AIC
#   3. For 10 new points, compute adjusted and unadjusted intervals
#   4. Check coverage (whether true value falls in interval)
# Returns: list with selected model, coverage indicators, true y values
one_rep <- function(j) {
  
  # ----- Simulation parameters -----
  p <- 8                    # Number of predictors
  B <- rep(0, p)
  B[1:3] <- 1:3             # True coefficients: only first 3 are nonzero
  statistic  <- "aic"       # Model selection criterion: "aic", "aicc", "bic"
  alpha      <- 0.05        # Significance level for intervals
  sigmaKnown <- FALSE       # Treat variance as unknown (estimate from data)
  
  cor_matrix <- diag(p)     # Independent design
  n <- 50                   # Sample size
  
  # ----- Generate 10 new points for prediction (fixed across replications) -----
  set.seed(1)
  new_X <- MASS::mvrnorm(10, mu = rep(0, p), Sigma = cor_matrix)
  
  # ----- Enumerate all possible competing models -----
  # (Binary matrix: each row is a candidate model)
  z <- expand.grid(rep(list(0:1), p))[-1, ]
  all_compete_models <- (matrix(unlist(z), ncol = p, byrow = FALSE)) == 1
  
  # ----- Generate training data for this replication -----
  set.seed(j)
  
  X  <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
  df <- data.frame(X)
  y  <- X %*% B + rnorm(nrow(X), 0, 1)  # Error variance = 1
  df$Y <- y
  
  all_vars <- paste0("X", 1:p)
  
  # ----- Best-subset model selection -----
  models <- leaps::regsubsets(
    Y ~ .,
    data = df,
    nvmax = p,
    method = "exhaustive",
    nbest = 1,
    intercept = TRUE
  )
  reg_summary <- summary(models)
  
  # choose subset size by AIC/AICC/BIC
  if (statistic == "aic") {
    aics <- rep(NA_real_, ncol(X))
    for (i in 1:ncol(X)) {
      vars <- all_vars[reg_summary$which[i, -1]]
      fmla <- as.formula(paste("Y ~", paste(vars, collapse = "+")))
      aics[i] <- AIC(lm(fmla, data = df))
    }
    size <- which.min(aics)
    
  } else if (statistic == "bic") {
    size <- which.min(reg_summary$bic)
    
  } else if (statistic == "aicc") {
    aiccs <- rep(NA_real_, ncol(X))
    for (i in 1:ncol(X)) {
      vars <- all_vars[reg_summary$which[i, -1]]
      fmla <- as.formula(paste("Y ~", paste(vars, collapse = "+")))
      k <- length(vars)
      aiccs[i] <- AIC(lm(fmla, data = df)) + 2 * k * (k + 1) / (n - k - 1)
    }
    size <- which.min(aiccs)
    
  } else stop("Invalid statistic!")
  
  selected <- reg_summary$which[size, -1]
  
  # ----- Compute coverage for 10 new points -----
  # Preallocate storage for speed
  y_list  <- numeric(10)  # True y values at new points
  c1_list <- numeric(10)  # Coverage indicator for adjusted intervals
  c2_list <- numeric(10)  # Coverage indicator for naive intervals
  
  for (i in 1:10) {
    new_x <- new_X[i, ]
    newy  <- as.numeric(new_x %*% B)
    y_list[i] <- newy
    
    interval1 <- saturated_interval(
      X_dt = X,
      y = y,
      selected = selected,
      contrast = "new",
      all_compete_models = all_compete_models,
      new_xpoint = new_x,
      statistic = statistic,
      true_sigma = NA,
      alpha = alpha
    )
    c1_list[i] <- as.numeric(newy > interval1[1] && newy < interval1[2])
    
    interval2 <- naive_interval(
      X_dt = X,
      y = y,
      selected = selected,
      contrast = "new",
      all_compete_models = all_compete_models,
      new_xpoint = new_x,
      statistic = statistic,
      true_sigma = NA,
      alpha = alpha
    )
    c2_list[i] <- as.numeric(newy > interval2[1] && newy < interval2[2])
  }
  
  list(
    selected = selected,
    c1 = c1_list,
    c2 = c2_list,
    ynew = y_list
  )
}


## --- Test single replication ---
one_rep(1)

## --- Run 500 replications in parallel ---
res <- mclapply(1:500, one_rep, mc.cores = max(1, detectCores() - 1))

## --- Extract results from all replications ---
se_List   <- lapply(res, `[[`, "selected")
c1_List   <- lapply(res, `[[`, "c1")
c2_List   <- lapply(res, `[[`, "c2")
ynew_List <- lapply(res, `[[`, "ynew")

## --- Compute empirical coverage rates ---
# Adjusted (selective) intervals
mat <- do.call(rbind, c1_List)  # dim: 500 x 10
adjusted = colSums(mat)/500     # Coverage for each of 10 new points


# Unadjusted (naive) intervals
mat <- do.call(rbind, c2_List)  # dim: 500 x 10
unadjusted = colSums(mat)/500   # Coverage for each of 10 new points

## --- Prepare data for plotting ---
d1 = data.frame(adjusted)
d1$method = 'adjusted'
d1$x = 1:10
colnames(d1) = c('Coverage','Method','x')


d2 = data.frame(unadjusted)
d2$method = 'unadjusted'
d2$x = 1:10
colnames(d2) = c('Coverage','Method','x')

result = rbind(d1,d2)

## --- Create coverage plot with confidence bands ---
# Prepare data frame
df = result
# df has columns: Coverage, Method, x
# Ensure proper data types
df$x <- factor(df$x)           # Treat x as discrete (1..10) for point plots
df$Method <- factor(df$Method)

# Add 95% confidence bands (binomial proportion, normal approximation)
df$upper = df$Coverage + sqrt(df$Coverage*(1-df$Coverage)/500)*1.99
df$lower = df$Coverage - sqrt(df$Coverage*(1-df$Coverage)/500)*1.99






ggplot(df, aes(x = x, y = Coverage, color = Method)) +
  geom_errorbar(
    aes(ymin = lower,
        ymax = upper),
    width = 0.15,                        # width of bar ends
    position = position_dodge(width = 0.4)
  ) +
  geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  #geom_line(aes(group = Method), position = position_dodge(width = 0.3)) +
  labs(x = "x", y = "Coverage") +
  theme(legend.position = "right")

## --- Stratified coverage analysis by selected model size ---
# Examine coverage conditional on how many predictors were selected
m_size = sapply(se_List, sum)
table(m_size)  # Frequency of each model size

mean(unlist(c1_List[m_size == 3]))
mean(unlist(c2_List[m_size == 3]))

mean(unlist(c1_List[m_size == 4]))
mean(unlist(c2_List[m_size == 4]))

mean(unlist(c1_List[m_size == 5]))
mean(unlist(c2_List[m_size == 5]))

mean(unlist(c1_List[m_size == 6]))
mean(unlist(c2_List[m_size == 6]))

