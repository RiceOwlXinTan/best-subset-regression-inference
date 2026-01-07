

################################################################################
# Application example: post-selection inference on `fpp3::us_change`
#
# This script demonstrates how to use the functions in `source code.R` to:
#   1) Select a linear regression model via an information criterion (AIC/BIC/AICc)
#   2) Compute p-values under different inference approaches:
#        - saturated_p(): selective p-value under the saturated model
#        - naive_p(): classical p-value ignoring selection
#        - selective_P_classic_fast(): Monte Carlo selective p-value (re-selection)
#   3) Compute confidence intervals via inversion:
#        - saturated_interval(), naive_interval()
#
# Data
# - Uses `us_change` from the `fpp3` package.
# - Predictors: columns 3:6
# - Response: `Consumption`
#
# Important notes / assumptions
# - This file sources `source code.R`, which currently contains a sanity-check
#   block that runs on source (plots + quick checks). If you only want function
#   definitions, consider commenting out that block.
# - The variables `sigmaKnown` and `new_x` are referenced below but not defined
#   in this script. In the loops here we test coefficients (contrast1 = 1:4),
#   so `new_x` is not used unless you set contrast1 = 'new'. You can set:
#     sigmaKnown <- TRUE
#     new_x <- rep(NA_real_, 4)
#   before running if you want a clean, self-contained example.
# - In the selective p-value call, the argument name should match the function
#   signature (`draws`, not `draw`). This script keeps the original line to
#   match the submission code; adjust locally if you run into an error.
################################################################################

source('source code.R')
library(fpp3)
library(intervals)
library(intervals) # duplicated in original script; safe but redundant
library(nleqslv)
library(leaps)
library(natural)
library(MASS)
library(parallel)
library(foreach)
library(doParallel)
library(bigstatsr)


## --- Step 1: set up design matrix and selection criterion ---
n = 198
statistic <- "aic"
X = us_change[,3:6]
df = data.frame(X)
df$Y= us_change$Consumption
p = 4
## --- Step 2: best-subset model selection via information criterion ---
all_vars <- c(paste("X", 1:p, sep=""))
colnames(df)[1:4] = all_vars


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



## --- Selected model and (optional) sigma estimate ---
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



## --- Step 3: p-values ---
## Saturated-model selective p-values (one per coefficient in the selected set)

p_l = c()
for (j in 1:4){
  p = saturated_p(X_dt = as.matrix(us_change[,3:6]), 
                   y = as.matrix(us_change$Consumption), 
                   selected,
                   contrast1 = j,
                   contrast2 = NA,
                   alls = all_compete_models, 
                   new_xpoint = new_x, 
                   statistic='aic',
                   true_sigma ='full',
                   null_value = 0)
  p_l = c(p_l,p)
  
}
round(p_l, digits = 4)


## Naive p-values (ignoring selection)
p_l = c()
for (j in 1:4){
  p = naive_p(X_dt = as.matrix(us_change[,3:6]),
          y = as.matrix(us_change$Consumption),
          selected,
          contrast1 = j,
          contrast2 = NA,
          new_xpoint = new_x,
          statistic='aic',
          null_value = 0,
          true_sigma = NA,
          make_plot = TRUE)
  p_l = c(p_l,p)
  
}
round(p_l, digits = 4)


## Selective p-values via Monte Carlo re-selection (classic-fast)
p_l = c()
for (j in 1:4){
  p = selective_P_classic_fast(X_dt = as.matrix(us_change[,3:6]),
                  y = as.matrix(us_change$Consumption),
                  selected,
                  contrast1 = j,
                  contrast2 = NA,
                  alls = all_compete_models,
                  new_xpoint = new_x,
                  statistic='aic',
                  draw = 1000,
                  true_sigma =NA,
                  null_value = 0)
  p
  p_l = c(p_l,p)

}
round(p_l, digits = 4)


#0.000 0.024 0.000 0.104


## --- Step 4: confidence intervals ---
## Saturated intervals (invert saturated_p)
for (j in 1:4){
  interval = saturated_interval(X_dt = as.matrix(us_change[,3:6]),
                     y = as.matrix(us_change$Consumption),
                     selected,
                     contrast = j,
                     all_compete_models = all_compete_models,
                     new_x = NA,
                     statistic='aic',
                     true_sigma = NA,
                     alpha = 0.05)
  print( interval)
}


## Naive intervals (ignore selection)
for (j in 1:4){
  interval = naive_interval(X_dt = as.matrix(us_change[,3:6]),
                 y = as.matrix(us_change$Consumption),
                 selected,
                 contrast = j,
                 all_compete_models = all_compete_models,
                 new_x = NA,
                 statistic='aic',
                 true_sigma = NA,
                 alpha = 0.05)
  print( interval)
}




