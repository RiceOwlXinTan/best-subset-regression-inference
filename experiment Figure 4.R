
source('source code.R')
#load packages
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
library(tidyr)
library(ggplot2)
library(dplyr)
library(patchwork)


selective_exp = function(b2,
                          n_samples = 50,
                          noise_level = 1,
                          noise = 'known',
                          p_samples = 100){
  
  p_list1 = c()
  p_list2 = c()
  p_list3 = c()
  s_list  = c()
  
  n_selected = 0
  iteration  = 0
  
  
  while (n_selected<p_samples) {
    
    iteration = iteration + 1
    
    p = 5 # num of columns in X
    B = rep(0,p)
    B[1]    = 1 
    B[2]    = b2
    B[3]    = 2
    B[4]    = 0
    B[5]    = 0.5
    
    n = n_samples
    statistic <- "aic" #must be one of c("aic", "aicc", "bic")
    alpha = 0.05
    sigmaKnown=F#assume noise level sigma is unknown
    
    #===Data generation==============
    # correlation structure of X: AR(1) with rho = 0.5
    cor_matrix=ar1_matrix(p, 0.5)
    
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
    z=expand.grid(rep(list(0:1),p))[-1,]
    all_compete_models=(matrix(unlist(z), ncol = p, byrow = F))==1
    selected = reg_summary$which[size,-1]
    
    if (selected[2] == TRUE){
      
      if (noise == 'known'){
        sigma_value = noise_level
      } else{
        sigma_value = NA
      }
      
      s_list = rbind(s_list,selected)
      
      n_selected =  n_selected + 1
      
      p = saturated_p(X_dt = X, 
                      y, 
                      selected,
                      contrast1 = sum(selected[1:2]),
                      contrast2 = NA,
                      alls = all_compete_models, 
                      new_xpoint = new_x, 
                      statistic='aic',
                      true_sigma = sigma_value,
                      null_value = 0)
      p_list1 = c(p_list1,p)
      
      
      p = naive_p(X_dt = X,
                  y, 
                  selected,
                  contrast1 = sum(selected[1:2]),
                  contrast2 = NA,
                  new_xpoint = new_x, 
                  statistic='aic', 
                  true_sigma = sigma_value,
                  null_value = 0)
      p_list2 = c(p_list2,p)
      
      
      p = selective_P_classic_fast(X_dt = X,
                                   y,
                                   selected,
                                   contrast1 = sum(selected[1:2]),
                                   contrast2 = NA,
                                   alls = all_compete_models,
                                   new_xpoint = new_x,
                                   statistic='aic',
                                   true_sigma = sigma_value,
                                   draws = 500,
                                   null_value = 0)
      p_list3 = c(p_list3,p)
      
      
    }
    if (iteration > 1000000){
      break
    }
  }
  return (list(p_list1,
               p_list2, 
               p_list3))
}



b_candidates <-  seq(-1.5,1.5,0.1)

# how many cores to use:
ncores <- detectCores() - 2

# run in parallel:
res_list1 <- mclapply(b_candidates, function(b) {
  out <- selective_exp (b,
                        n_sample = 30,
                        p_samples = 1000,
                        noise = 'unknown',
                        noise_level = 1)
  c(
    rejection1 = mean(out[[1]] < 0.05),
    rejection2 = mean(out[[2]] < 0.05),
    rejection3 = mean(out[[3]] < 0.05)
    )
  }, 
  mc.cores = ncores)



res_list1


# load precomputed result
load('my_vars_new.RData')






########################################Figure 4 known variance

# combine results into two vectors
res_mat   <- do.call(rbind, res_list1)


rejection1 <- res_mat[, "rejection1"]
rejection2 <- res_mat[, "rejection2"]
rejection3 <- res_mat[, "rejection3"]


# base plot
plot(b_candidates,  rejection1, type = "l", xlab = "b", ylab = "Rejection rate", ylim = c(0,1))
lines(b_candidates, rejection2, col = "red")
lines(b_candidates, rejection3, col = "blue")
abline(h = 0.05)
legend("bottomright", c("saturated","naive"), col = c("black","red"), lty = 1)



df = data.frame(b_candidates, rejection1, rejection2, rejection3)
colnames(df) = c('b','saturated', 'naive', 'selective')

# Long format
df_long <- pivot_longer(
  df,
  cols = c(saturated, naive, selective),
  names_to = "method",
  values_to = "rejection"
)

# Single labeling variable (4 combos)
df_long$label <- factor(
  dplyr::recode(df_long$method,
                saturated = "Saturated",
                naive = "Naive",
                selective    = "Selective"
  ),
  levels = c("Saturated",
             "Naive",
             "Selective")
)



# Plot: one legend with 4 entries 
p1 = ggplot(df_long, aes(x = b, y = rejection, color = label, linetype = label)) +
  geom_line(size = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed",color = 'blue') +
  scale_color_manual(values = c(
    "Saturated"  = "black",
    "Naive" = "blue",
    "Selective"     = "red"
  )) +
  scale_linetype_manual(values = c(
    "Saturated"  = "solid",
    "Naive" = "solid",
    "Selective"     = "solid"
  )) +
  labs(x =expression(beta[2]) , y = "Rejection rate", color = NULL, linetype = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  #theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")+
  ggtitle('Variance Known')


# add confidence interval for the rejection
df_ci <- df_long %>%
  mutate(
    lower = pmax(0, rejection - 1.98*sqrt(rejection*(1-rejection)/500)),
    upper = pmin(1, rejection + 1.98*sqrt(rejection*(1-rejection)/500))
  )


p1 <- ggplot() +
  # confidence band
  geom_ribbon(
    data = df_ci,
    aes(x = b, ymin = lower, ymax = upper, fill = label),
    alpha = 0.2,            # transparency
    inherit.aes = FALSE
  ) +
  # original lines
  geom_line(
    data = df_long,
    aes(x = b, y = rejection, color = label, linetype = label),
    size = 0.5
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +
  scale_color_manual(values = c(
    "Saturated" = "black",
    "Naive"     = "blue",
    "Selective" = "red"
  )) +
  scale_fill_manual(values = c(   # match band colors to lines
    "Saturated" = "black",
    "Naive"     = "blue",
    "Selective" = "red"
  )) +
  scale_linetype_manual(values = c(
    "Saturated" = "solid",
    "Naive"     = "solid",
    "Selective" = "solid"
  )) +
  labs(
    x = expression(beta[2]),
    y = "Rejection rate",
    color = NULL,
    linetype = NULL,
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = "bottom") +
  ggtitle("Variance Known")


######################################################## Figure 4 unknown variance


# combine results into two vectors
res_mat   <- do.call(rbind, res_list2)


rejection1 <- res_mat[, "rejection1"]
rejection2 <- res_mat[, "rejection2"]
rejection3 <- res_mat[, "rejection3"]

# plot
plot(b_candidates,  rejection1, type = "l", xlab = "b", ylab = "Rejection rate", ylim = c(0,1))
lines(b_candidates, rejection2, col = "red")
lines(b_candidates, rejection3, col = "blue")
abline(h = 0.05)
legend("bottomright", c("saturated","naive"), col = c("black","red"), lty = 1)


df = data.frame(b_candidates, rejection1, rejection2, rejection3)
colnames(df) = c('b','saturated', 'naive', 'selective')

# Long format
df_long <- pivot_longer(
  df,
  cols = c(saturated, naive, selective),
  names_to = "method",
  values_to = "rejection"
)

# Single labeling variable (4 combos)
df_long$label <- factor(
  dplyr::recode(df_long$method,
                saturated = "Saturated",
                naive = "Naive",
                selective    = "Selective"
  ),
  levels = c("Saturated",
             "Naive",
             "Selective")
)

# Plot: one legend with 4 entries
p2 = ggplot(df_long, aes(x = b, y = rejection, color = label, linetype = label)) +
  geom_line(size = 0.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed",color = 'blue') +
  scale_color_manual(values = c(
    "Saturated"  = "black",
    "Naive" = "blue",
    "Selective"     = "red"
  )) +
  scale_linetype_manual(values = c(
    "Saturated"  = "solid",
    "Naive" = "solid",
    "Selective"     = "solid"
  )) +
  labs(x =expression(beta[2]) , y = "Rejection rate", color = NULL, linetype = NULL) +
  coord_cartesian(ylim = c(0, 1)) +
  #theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")+
  ggtitle('Variance Unknown')

df_ci <- df_long %>%
  mutate(
    lower = pmax(0, rejection - 1.98*sqrt(rejection*(1-rejection)/500)),
    upper = pmin(1, rejection + 1.98*sqrt(rejection*(1-rejection)/500))
  )

p2 <- ggplot() +
  # confidence band
  geom_ribbon(
    data = df_ci,
    aes(x = b, ymin = lower, ymax = upper, fill = label),
    alpha = 0.2,            # transparency
    inherit.aes = FALSE
  ) +
  # original lines
  geom_line(
    data = df_long,
    aes(x = b, y = rejection, color = label, linetype = label),
    size = 0.5
  ) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "blue") +
  scale_color_manual(values = c(
    "Saturated" = "black",
    "Naive"     = "blue",
    "Selective" = "red"
  )) +
  scale_fill_manual(values = c(   # match band colors to lines
    "Saturated" = "black",
    "Naive"     = "blue",
    "Selective" = "red"
  )) +
  scale_linetype_manual(values = c(
    "Saturated" = "solid",
    "Naive"     = "solid",
    "Selective" = "solid"
  )) +
  labs(
    x = expression(beta[2]),
    y = "Rejection rate",
    color = NULL,
    linetype = NULL,
    fill = NULL
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme(legend.position = "bottom") +
  ggtitle("Variance Unknown")


#################################### Final figure

(p1 | p2) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")


