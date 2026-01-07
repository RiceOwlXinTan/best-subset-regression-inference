################################################################################
# Experiment: Figures 2, 3, and 5 from the paper
#
# Purpose:
#   Generates three key figures illustrating post-selection inference:
#
#   Figure 2: Null distributions under different inference methods
#     - Shows how saturated, naive, and selective approaches differ
#     - Visualizes density/mass functions for selected coefficients
#
#   Figure 3: Type I error control across effect sizes
#     - Rejection rates as a function of beta_2
#     - Compares known vs unknown variance settings
#
#   Figure 5: Conditional type I error by model selected
#     - Stratifies rejection rates by which model was chosen
#     - Demonstrates importance of conditioning on selection event
#
# Setup:
#   - p = 10 predictors (true model has first 3 active)
#   - Design: independent (Figure 2/3) or AR(1) with rho=0.5 (Figure 5)
#   - Model selection via AIC
#
# Note: Heavy computation - uses parallel processing via mclapply().
#       Precomputed results can be loaded from .rds files.
################################################################################

rm(list = ls())
source('source code.R')

## --- Load required packages ---
library(intervals)
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
library(parallel)
library(foreach)
library(doParallel)
library(foreach)
library(doParallel)
library(tidyr)
library(ggplot2)
library(parallel)
library(latex2exp)
library(patchwork)
library(fpp3)
library(final)

## ============================================================================
##  SECTION 1: Single example data for Figure 2 null distributions
## ============================================================================

set.seed(4)

## --- Setup for single example dataset ---
p=10 # Number of predictors
B=rep(0,p)
B[1:3]=1:3 #ground truth model: only the first 3 predictors have nonzero effect on y

statistic <- "aic" #must be one of c("aic", "aicc", "bic")
alpha=0.05
sigmaKnown=F#assume noise level sigma is unknown


#===Data generation==============
# correlation structure of X: AR(1) with rho = 0.5
cor_matrix=diag(p)

#generate data used in model fitting

n = 50
X=MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
df=data.frame(X)
y=X%*%B + rnorm(nrow(X), 0, 1)
df$Y=y
#generate a new data point
new_x=MASS::mvrnorm(1, mu = rep(0, p), Sigma = cor_matrix)


#===Go through the information-based model selection procedure========
all_vars <- c(paste("X", 1:p, sep=""))
names(new_x)=all_vars

models <- leaps::regsubsets(Y~., data = df, nvmax = p, method = 'exhaustive', nbest = 1, intercept=TRUE)
reg_summary = summary(models)
if (statistic == 'aic') {
  aics <- rep(NA,ncol(X))
  for (i in 1:ncol(X)) {
    vars <- all_vars[reg_summary$which[i,-1]]
    # Fit the linear regression with the variables in subset.
    fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
    est_fit <- lm(fmla, data = df)
    aics[i] <- AIC(est_fit)
  }
  size = which.min(aics)
} else if (statistic == 'bic') {
  size = which.min(reg_summary$bic)
} else if(statistic=='aicc'){
  aiccs <- rep(NA,ncol(X))
  for (i in 1:ncol(X)) {
    vars <- all_vars[reg_summary$which[i,-1]]
    # Fit the linear regression with the variables in subset.
    fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
    est_fit <- lm(fmla, data = df)
    aiccs[i] <- AIC(est_fit)+2*length(vars)*(length(vars)+1)/(n-length(vars)-1)
  }
  size = which.min(aiccs)
}else {
  print("Invalid statistic!")
  break;
}
#selected model
vars <- all_vars[reg_summary$which[size,-1]]#selected variables
fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
est_fit <- lm(fmla, data = df)
#estimated noise level, if unknown
if(!sigmaKnown){
  sigmahat=summary(est_fit)$sigma
}

#===Get post-selection corrected CI for predicted mean at new_x
#all possible models (use 1 to indicate a variable is selected and 0 otherwise)
z=expand.grid(rep(list(0:1),p))[-1,]
all_compete_models=(matrix(unlist(z), ncol = p, byrow = F))==1
selected = reg_summary$which[size,-1]

## ============================================================================
##  SECTION 2: Generate Figure 2 - Null distributions
## ============================================================================

## --- Test single p-value computation ---
selective_P_classic_fast(X_dt = X,
                         y,
                         selected,
                         contrast1 = 5,
                         contrast2 = NA,
                         alls = all_compete_models,
                         new_xpoint = new_x,
                         statistic='aic',
                         draw = 100,
                         true_sigma = NA,
                         null_value = 0)





## --- Compute null distributions for selected coefficients ---
# Grid: test 3 coefficients (1, 4, 5) x 2 variance settings (known/unknown)
grid <- expand.grid(i = c(1,4,5), s = c(NA,1), KEEP.OUT.ATTRS = FALSE)

res <- mclapply(seq_len(nrow(grid)), function(k) {
  i <- grid$i[k]; s <- grid$s[k]
  null <- if (i == 1) 1 else 0
  ap   <- if (i == 1) 1 else if (i == 4) 4 else 9
  
  p1 <- saturated_p(X_dt = X, y, selected, contrast1 = i, contrast2 = NA,
                    alls = all_compete_models, new_xpoint = new_x,
                    statistic="aic", null_value=null, true_sigma=s, dist_null=TRUE)
  
  p2 <- naive_p(X_dt = X, y, selected, contrast1 = i, contrast2 = NA,
                new_xpoint = new_x, statistic="aic", null_value=null,
                true_sigma=s, make_plot=TRUE, dist_null=TRUE)
  
  p3 <- selective_P_classic_fast(X_dt = X, y, selected, contrast1 = i, contrast2 = NA,
                                 actual_p=ap, alls=all_compete_models, new_xpoint=new_x,
                                 statistic="aic", draw=1000, true_sigma=s,
                                 null_value=null, dist_null=TRUE)
  
  list(p1=p1, p2=p2, p3=p3)
}, mc.cores = detectCores())


## --- Save/load precomputed results ---
#saveRDS(res,'pivot.rds')
res = readRDS('pivot.rds')

## --- Extract density curves for plotting ---
sat_l <- lapply(res, `[[`, "p1")
nai_l <- lapply(res, `[[`, "p2")
sel_l <- lapply(res, `[[`, "p3")


x_l  <- seq(-3, 3, length.out = 1000)


s_d = density(sel_l[[2]])
x <- s_d$x
y <- s_d$y
x_new <- seq(min(x), max(x), by = 0.01001001)
y_new <- approx(x, y, xout = x_new)$y

df <- rbind(
  data.frame(x = x_l,  y = sat_l[[2]], method = "Saturated"),
  data.frame(x = x_l,  y = nai_l[[2]], method = "Naive"),
  data.frame(x = x_new, y = y_new * 1.6, method = "Selective")
)

p1 = ggplot(df, aes(x = x, y = y, fill = method, color = method)) +
  geom_ribbon(aes(ymin = 0, ymax = y),
              alpha = 0.35,
              color = NA) +
  xlim(-1,1)+
  ylim(0,14)+
  geom_line(linewidth = 0.3) +
  scale_color_manual(values = c(
    "Saturated"    = "black",
    "Naive"        = "blue",
    "Selective" = "red"
  )) +
  scale_fill_manual(values = c(
    "Saturated"    = "black",
    "Naive"        = "blue",
    "Selective" = "red"
  ))+
  ylab('Density')+
  xlab(TeX("Null distribution for $\\beta_4$"))



s_d = density(sel_l[[3]])
x <- s_d$x
y <- s_d$y
x_new <- seq(min(x), max(x), by = 0.01001001)
y_new <- approx(x, y, xout = x_new)$y

df <- rbind(
  data.frame(x = x_l,  y = sat_l[[3]], method = "Saturated"),
  data.frame(x = x_l,  y = nai_l[[3]], method = "Naive"),
  data.frame(x = x_new, y = y_new * 1.6, method = "Selective")
)

p2 = ggplot(df, aes(x = x, y = y, fill = method, color = method)) +
  geom_ribbon(aes(ymin = 0, ymax = y),
              alpha = 0.35,
              color = NA) +
  geom_line(linewidth = 0.3) +
  xlim(-1,1)+
  ylim(0,14)+
  scale_color_manual(values = c(
    "Saturated"    = "black",
    "Naive"        = "blue",
    "Selective" = "red"
  )) +
  scale_fill_manual(values = c(
    "Saturated"    = "black",
    "Naive"        = "blue",
    "Selective" = "red"
  ))+
  ylab('Density')+
  xlab(TeX("Null distribution for $\\beta_9$"))

p2

p1

(p1 | p2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")


## ============================================================================
##  SECTION 3: Generate Figure 3 - Type I error across effect sizes
## ============================================================================

## --- Main simulation function for Figure 3 ---
# Generates data, selects model, computes p-values under different methods
saturated_exp1 = function(b2,
                          n_samples = 50,
                          noise_level = 1,
                          p_samples = 100){
  
  p_list1 = c()
  p_list2 = c()
  p_list3 = c()
  p_list4 = c()
  s_list  = c()
  n_selected = 0
  iteration  = 0
  
  
  while (n_selected<p_samples) {
    
    iteration = iteration + 1
    
    
    p = 10 # num of columns in X
    B = rep(0,p)
    B[1]    = 1
    B[2]    = 2
    B[3]    = 3
    B[4]    = 0
    B[5]    = 0
    B[6]    = 0
    B[7]    = 0
    B[8]    = 0
    B[9]    = 0
    B[10]   = 0
    
    
    n = n_samples
    statistic <- "aic" #must be one of c("aic", "aicc", "bic")
    alpha = 0.05
    sigmaKnown=F#assume noise level sigma is unknown
    
    #===Data generation==============
    # correlation structure of X: AR(1) with rho = 0.5
    cor_matrix= diag(p)*0.5
      # ar1_matrix(p, 0.5)
    
    #generate data used in model fitting
    set.seed(iteration)
    
    X=MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = cor_matrix)
    df=data.frame(X)
    y=X%*%B + rnorm(nrow(X), 0, noise_level)
    df$Y=y
    #generate a new data point
    new_x=MASS::mvrnorm(1, mu = rep(0, p), Sigma = cor_matrix)
    
    #===Go through the information-based model selection procedure========
    all_vars <- c(paste("X", 1:p, sep=""))
    names(new_x)=all_vars
    
    models <- leaps::regsubsets(Y~., 
                                data = df, 
                                nvmax = p, 
                                method = 'exhaustive', 
                                nbest = 1, 
                                intercept=TRUE)
    reg_summary = summary(models)
    if (statistic == 'aic') {
      aics <- rep(NA,ncol(X))
      for (i in 1:ncol(X)) {
        vars <- all_vars[reg_summary$which[i,-1]]
        # Fit the linear regression with the variables in subset.
        fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
        est_fit <- lm(fmla, data = df)
        aics[i] <- AIC(est_fit)
        }
      size = which.min(aics)
      } 
    else if (statistic == 'bic') {
        size = which.min(reg_summary$bic)
        } 
    else if (statistic=='aicc') {
          aiccs <- rep(NA,ncol(X))
          for (i in 1:ncol(X)) {
            vars <- all_vars[reg_summary$which[i,-1]]
            # Fit the linear regression with the variables in subset.
            fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
            est_fit <- lm(fmla, data = df)
            aiccs[i] <- AIC(est_fit)+2*length(vars)*(length(vars)+1)/(n-length(vars)-1)
              }
          size = which.min(aiccs)
          }
    else {
      print("Invalid statistic!")
      break;
    }
    
    #selected model
    vars <- all_vars[reg_summary$which[size,-1]]#selected variables
    fmla <- as.formula(paste("Y~", paste(vars, collapse="+")))
    est_fit <- lm(fmla, data = df)
    #estimated noise level, if unknown
    if(!sigmaKnown){
      sigmahat=summary(est_fit)$sigma
    }
    
    #===Get post-selection corrected CI for predicted mean at new_x
    #all possible models (use 1 to indicate a variable is selected and 0 otherwise)
    all_compete_models=(matrix(unlist(z), ncol = p, byrow = F))==1
    selected = reg_summary$which[size,-1]
    
    # match = which(apply(all_compete_models, 1, function(r) all(r == selected)))
    # all_compete_models = all_compete_models[-match,]
    
    
    if (selected[2] == TRUE){
      
      s_list = rbind(s_list,selected)
      
      n_selected =  n_selected + 1
      
      p_val = saturated_p(X_dt = X,
                      y,
                      selected,
                      contrast1 = sum(selected[1:2]),
                      contrast2 = NA,
                      alls = all_compete_models,
                      new_xpoint = new_x,
                      statistic='aic',
                      true_sigma = noise_level,
                      null_value = 0)
      p_list1 = c(p_list1,p_val)
      
      p_val = naive_p(X_dt = X,
                      y, 
                      selected,
                      contrast1 = sum(selected[1:2]),
                      contrast2 = NA,
                      new_xpoint = new_x, 
                      statistic='aic', 
                      true_sigma = noise_level,
                      null_value = 0)
      p_list2 = c(p_list2,p_val)
      
      
      # p = selective_P_classic_fast(X_dt = X,
      #                 y,
      #                 selected,
      #                 contrast1 = sum(selected[1:2]),
      #                 contrast2 = NA,
      #                 alls = all_compete_models,
      #                 new_xpoint = new_x,
      #                 statistic='aic',
      #                 true_sigma = noise_level,
      #                 draws = 500,
      #                 null_value = 0)
     
      p_val = saturated_p(X_dt = X,
                      y,
                      selected,
                      contrast1 = sum(selected[1:2]),
                      contrast2 = NA,
                      alls = all_compete_models,
                      new_xpoint = new_x,
                      statistic='aic',
                      true_sigma = NA,
                      null_value = 0)
      p_list3 = c(p_list3,p_val)
      
      
      p_val = naive_p(X_dt = X,
                  y,
                  selected,
                  contrast1 = sum(selected[1:2]),
                  contrast2 = NA,
                  new_xpoint = new_x,
                  statistic='aic',
                  true_sigma = NA,
                  null_value = 0)
      p_list4 = c(p_list4, p_val)
      
    }
    if (iteration > 1000000){
      break
    }
  }
  return (list(p_list1,
               p_list2, 
               p_list3,
               p_list4))
}

## --- Test single run at null (beta_2 = 0) ---
result = saturated_exp1(b = 0,n_sample = 50, p_samples = 100,  noise_level = 1)

## --- Run simulations across range of beta_2 values ---
b_candidates <-  seq(-1.5,1.5,0.1)


all(all_compete_models[s,]-selected>=0)

s = 3
# how many cores to use:
ncores <- detectCores() - 2

# run in parallel:
res_list <- mclapply(b_candidates, function(b) {
  out <- saturated_exp1(b,n_sample = 100, p_samples = 1000,  noise_level = 1)
  c(
    rejection1 = mean(out[[1]] < 0.05),
    rejection2 = mean(out[[2]] < 0.05),
    rejection3 = mean(out[[3]] < 0.05),
    rejection4 = mean(out[[4]] < 0.05)
    )
  }, 
  mc.cores = ncores)


# combine results into two vectors
res_mat   <- do.call(rbind, res_list)



rejection1 <- res_mat[, "rejection1"]
rejection2 <- res_mat[, "rejection2"]
rejection3 <- res_mat[, "rejection3"]
rejection4 <- res_mat[, "rejection4"]
# plot
plot(b_candidates,  rejection1, type = "l", xlab = "b", ylab = "Rejection rate", ylim = c(0,1))
lines(b_candidates, rejection2, col = "red")
lines(b_candidates, rejection3, col = "black", lty = 2)
lines(b_candidates, rejection4, col = "red",   lty = 2)
abline(h = 0.05)
legend("bottomright", c("saturated","naive"), col = c("black","red"), lty = 1)

df = data.frame(b_candidates, rejection1, rejection2, rejection3, rejection4)
colnames(df) = c('b','adjusted1', 'naive1', 'adjusted2', 'naive2')

# Long format
df_long <- pivot_longer(
  df,
  cols = c(adjusted1, naive1, adjusted2, naive2),
  names_to = "method",
  values_to = "rejection"
)

# Single labeling variable (4 combos)
df_long$label <- factor(
  dplyr::recode(df_long$method,
                adjusted1 = "Adjusted (known σ)",
                adjusted2 = "Adjusted (unknown σ)",
                naive1    = "Naive (known σ)",
                naive2    = "Naive (unknown σ)"
  ),
  levels = c("Adjusted (known σ)",
             "Adjusted (unknown σ)",
             "Naive (known σ)", 
             "Naive (unknown σ)")
)

# Plot: one legend with 4 entries
ggplot(df_long, aes(x = b, y = rejection, color = label, linetype = label)) +
  geom_line(size = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed",color = 'blue') +
  scale_color_manual(values = c(
    "Adjusted (known σ)"  = "black",
    "Adjusted (unknown σ)"= "black",
    "Naive (known σ)"     = "red",
    "Naive (unknown σ)"   = "red"
  )) +
  scale_linetype_manual(values = c(
    "Adjusted (known σ)"  = "solid",
    "Adjusted (unknown σ)"= "dashed",
    "Naive (known σ)"     = "solid",
    "Naive (unknown σ)"   = "dashed"
  )) +
  labs(x =expression(beta[2]) , y = "Rejection rate", color = NULL, linetype = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  #theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")+
  title('n = 50')






## ============================================================================
##  SECTION 4: Generate Figure 5 - Conditional type I error by selected model
## ============================================================================

## --- Parallelized simulation function for Figure 5 ---
# Collects replications conditional on selecting a specific predictor
# More efficient than filtering post-hoc
saturated_exp4_par <- function(b2, 
                               test = 4,
                               n = 30, 
                               null_value = 0,
                               target_selected = 10,
                               chunk_size = 100, 
                               ncores = parallel::detectCores() - 1,
                               noise_level = 1,
                               statistic = "aic",
                               noise = 'known') {
  
  stopifnot(target_selected > 0, 
            chunk_size > 0, 
            statistic %in% c("aic","bic","aicc"))
  
  p        <- 10
  all_vars <- paste0("X", 1:p)
  
  s_rows  <- list()
  p_list  <- numeric(0)
  p_list2 <- numeric(0)
  p_list3 <- numeric(0)
  
  collected <- 0L
  iter0 <- 0L
 
  while (collected < target_selected) {
    iters <- (iter0 + 1L):(iter0 + chunk_size)
    
    res_chunk <- parallel::mclapply(iters, function(iter) {
      set.seed(iter)
      
      # ----- data gen -----
      
      B = rep(0,10)
      B[1]    = 1 
      B[2]    = 2
      B[3]    = 3
      B[4]    = 0
      # B[4]    = 0
      # B[5]    = 0.
      
      sigmaKnown=F#assume noise level sigma is unknown
      cor_matrix <- ar1_matrix(p, 0.5)                 # you must have this defined
      X <- MASS::mvrnorm(n = n, 
                         mu = rep(0, p), 
                         Sigma = cor_matrix)
      
      y=X%*%B + rnorm(nrow(X), 0, noise_level)
      
      df <- as.data.frame(X)
      names(df) <- all_vars
      df$Y <- y
      new_x <- setNames(as.numeric(MASS::mvrnorm(1, rep(0, p), cor_matrix)), all_vars)
      
      # ----- model selection -----
      models      <- leaps::regsubsets(Y ~ ., 
                                       data = df, 
                                       nvmax = p, 
                                       method = "exhaustive",
                                       nbest = 1, 
                                       intercept = TRUE)
      reg_summary <- summary(models)
      
      
      
      if (statistic == "aic") {
        aics <- vapply(1:p, function(i) {
          vars <- all_vars[reg_summary$which[i, -1]]
          fit  <- lm(as.formula(paste("Y~", paste(vars, collapse = "+"))), data = df)
          AIC(fit)
        }, numeric(1))
        size <- which.min(aics)
      } else if (statistic == "bic") {
        size <- which.min(reg_summary$bic)
      } else {
        aiccs <- vapply(1:p, function(i) {
          vars <- all_vars[reg_summary$which[i, -1]]
          fit  <- lm(as.formula(paste("Y~", paste(vars, collapse = "+"))), data = df)
          k <- length(vars) + 1
          AIC(fit) + 2*k*(k+1)/(n - k - 1)
        }, numeric(1))
        size <- which.min(aiccs)
      }
      

      
      selected <- reg_summary$which[size, -1]
      
      if (!isTRUE(selected[test])) return(NULL)
      
      if (noise == 'known'){
        sigma_value = noise_level
      } else{
        sigma_value = NA
      }
      
      
      # all competing models (logical matrix); avoids undefined 'z'
      alls <- reg_summary$which[, -1, drop = FALSE]
      
      # ----- your three p-values -----

      p_sat  <- saturated_p(X_dt = X,
                            y = y,
                            selected   = selected,
                            contrast1  = sum(selected[1:2]),
                            contrast2  = NA,
                            alls       = alls,
                            new_xpoint = new_x,
                            statistic  = statistic,
                            true_sigma = sigma_value,
                            null_value = 0)

      p_nv  <- naive_p(X_dt = X,
                       y    = y, 
                      selected   = selected,
                      contrast1  = sum(selected[1:test]), 
                      contrast2  = NA,
                      new_xpoint = new_x, 
                      statistic  = statistic, 
                      true_sigma = sigma_value,
                      null_value = null_value)
      
      p_fast <- selective_P_classic_fast(X_dt = X,
                                    y = y,
                                    selected = selected,
                                    contrast1 = sum(selected[1:2]),
                                    contrast2 = NA,
                                    actual_p = sum(selected[1:2]),
                                    alls = alls,
                                    new_xpoint = new_x,
                                    statistic = statistic,
                                    true_sigma = sigma_value,
                                    draws = 50,
                                    null_value = 0)
      
      list(selected = selected, 
           p_sat = p_sat ,
           p_naive = p_nv, 
           p_fast = p_fast )
    }, 
    mc.cores = ncores, 
    mc.preschedule = TRUE, 
    mc.set.seed = TRUE)
    
    # keep only successful selections
    res_chunk <- Filter(Negate(is.null), res_chunk)
    if (length(res_chunk)) {
      s_rows  <- c(s_rows,  lapply(res_chunk, `[[`, "selected"))
      p_list  <- c(p_list,  sapply(res_chunk, function(z) as.numeric(z$p_sat)[1]))
      p_list2 <- c(p_list2, sapply(res_chunk, function(z) as.numeric(z$p_naive)[1]))
      p_list3 <- c(p_list3, sapply(res_chunk, function(z) as.numeric(z$p_fast)[1]))
      collected <- length(p_list)
    }
    iter0 <- iter0 + chunk_size
  }
  
  s_mat <- do.call(rbind, s_rows)[seq_len(target_selected), , drop = FALSE]
  list(s_mat, 
       p_list[seq_len(target_selected)], 
       p_list2[seq_len(target_selected)], 
       p_list3[seq_len(target_selected)]
       )
}




## --- Run simulation (5000 replications conditional on X5 selected) ---
result_n = saturated_exp4_par(b2 = 0, 
                              n = 50, 
                              target_selected = 5000,
                              noise_level = 1, 
                              chunk_size = 1000,
                              null_value = 0,
                              test = 5,       # Condition on X5 being selected
                              noise = 'known')

result = result_n

## --- Save/load precomputed results ---
# saveRDS(result,'exp4.rds')
result = readRDS('exp4.rds')

## --- Examine distribution of selected models ---
combos <- apply(result[[1]], 1, paste0, collapse = "")
tab   <- sort(table(combos), decreasing = TRUE)
tab
indi = combos != names(tab[6])


par(mfrow = c(2, 2))



df <- data.frame(
  pval = result[[2]][indi]
)





## --- Compute rejection rates at different alpha levels ---
mean(result[[2]]<0.01)
mean(result[[3]]<0.05)
mean(result[[4]]<0.1)

## --- Construct Figure 5: rejection rates stratified by selected model ---
# Manually entered values from simulation results
# Format: c(condition, method, rejection_rate)

# Known variance results
r1 = c('all','naive',     0.2692)
r2 = c('all','saturated', 0.0604)
r3 = c('all','selected',  0.071)
r4 = c('m1', 'naive',     0.2473442)
r5 = c('m1', 'saturated', 0.05263992)
r6 = c('m1', 'selected',  0.05977383)
r7 = c('m2', 'naive',     0.2425604)
r8 = c('m2', 'saturated', 0.05053341)
r9 = c('m2', 'selected',  0.04891304)

u1 = c('all','naive',     0.287 )
u2 = c('all','saturated', 0.0913)
u3 = c('all','selected',  0.063)
u4 = c('m1', 'naive',     0.273021)
u5 = c('m1', 'saturated', 0.100271 )
u6 = c('m1', 'selected',  0.05815832)
u7 = c('m2', 'naive',     0.2771739 )
u8 = c('m2', 'saturated', 0.0932)
u9 = c('m2', 'selected',  0.05978261)

known   = as.data.frame(rbind(r1,r2,r3,r4,r5,r6,r7,r8,r9))
unknown = as.data.frame(rbind(u1,u2,u3,u4,u5,u6,u7,u8,u9))

colnames(known)   <- c("Model", "Method", "Rejection Rate")
colnames(unknown) <- c("Model", "Method", "Rejection Rate")

known$variance = 'known'
unknown$variance = 'unknown'
final = rbind(known, unknown)



# make sure things are factors (optional, but helps control order)
final$Model    <- factor(final$Model,   levels = c("all", "m1", "m2"))
final$Method   <- factor(final$Method,  levels = c("naive", "saturated", "selected"))
final$variance <- factor(final$variance, levels = c("known", "unknown"))
final$`Rejection Rate` = as.numeric(final$`Rejection Rate`)


final$upper = final$`Rejection Rate` + sqrt(final$`Rejection Rate`*(1-final$`Rejection Rate`)/500)*1.98
final$lower = final$`Rejection Rate` - sqrt(final$`Rejection Rate`*(1-final$`Rejection Rate`)/500)*1.98



ggplot(final,
       aes(x = Method,
           y = `Rejection Rate`,
           color = Model)) +
  geom_errorbar(
    aes(ymin = lower,
        ymax = upper),
    width = 0.15,                        # width of bar ends
    position = position_dodge(width = 0.4)
  ) +
  geom_point(position = position_dodge(width = 0.4), size = 3) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") + 
  facet_wrap(~ variance, nrow = 1) +
  scale_y_continuous(limits = c(0, 0.35)) +   
  labs(x = "Method",
       y = "Rejection rate",
       color = "Condition") 







