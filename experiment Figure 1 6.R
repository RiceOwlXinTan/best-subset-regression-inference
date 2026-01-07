# To replicate the simu results in sec 4 of the paper, with unknown sigma2.

#Overview: Three strategies for estimating sigma
#---strategy1: mse based on full model
#---strategy2: organic lasso
#---strategy3: mse based on aic-selected model

rm(list = ls())
source('source code.R')
#load packages
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

se_List = list()
c1_List = list()
c2_List = list()
ynew_List = list()

# ----------------------------
# Setup (outside the parallel loop)
# ----------------------------
# ----------------------------
# One replication as a function
# ----------------------------
one_rep <- function(j) {
  
  p <- 8
  B <- rep(0, p)
  B[1:3] <- 1:3
  statistic  <- "aic"   # "aic", "aicc", "bic"
  alpha      <- 0.05
  sigmaKnown <- FALSE
  
  cor_matrix <- diag(p)
  n <- 50
  
  set.seed(1)
  new_X <- MASS::mvrnorm(10, mu = rep(0, p), Sigma = cor_matrix)
  
  # all possible competing models depends only on p, compute once
  z <- expand.grid(rep(list(0:1), p))[-1, ]
  all_compete_models <- (matrix(unlist(z), ncol = p, byrow = FALSE)) == 1
  
  
  set.seed(j)
  
  X  <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
  df <- data.frame(X)
  y  <- X %*% B + rnorm(nrow(X), 0, 1)
  df$Y <- y
  
  all_vars <- paste0("X", 1:p)
  
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
  
  # preallocate for speed
  y_list  <- numeric(10)
  c1_list <- numeric(10)
  c2_list <- numeric(10)
  
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


one_rep(1)



res <- mclapply(1:500, one_rep, mc.cores = max(1, detectCores() - 1))

se_List   <- lapply(res, `[[`, "selected")
c1_List   <- lapply(res, `[[`, "c1")
c2_List   <- lapply(res, `[[`, "c2")
ynew_List <- lapply(res, `[[`, "ynew")

mat <- do.call(rbind, c1_List)  # dim: 100 x 10 (if 100 reps)
adjusted = colSums(mat)/500                    # length 10


mat <- do.call(rbind, c2_List)  # dim: 100 x 10 (if 100 reps)
unadjusted = colSums(mat)/500                    # length 10

d1 = data.frame(adjusted)
d1$method = 'adjusted'
d1$x = 1:10
colnames(d1) = c('Coverage','Method','x')


d2 = data.frame(unadjusted)
d2$method = 'unadjusted'
d2$x = 1:10
colnames(d2) = c('Coverage','Method','x')

result = rbind(d1,d2)




df = result
# df has columns: Coverage, Method, x
# make sure types are right
df$x <- factor(df$x)           # treat x as discrete 1..10 (common for point plots)
df$Method <- factor(df$Method)


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





m_size = sapply(se_List, sum)
table(m_size)

mean(unlist(c1_List[m_size == 3]))
mean(unlist(c2_List[m_size == 3]))

mean(unlist(c1_List[m_size == 4]))
mean(unlist(c2_List[m_size == 4]))

mean(unlist(c1_List[m_size == 5]))
mean(unlist(c2_List[m_size == 5]))

mean(unlist(c1_List[m_size == 6]))
mean(unlist(c2_List[m_size == 6]))

