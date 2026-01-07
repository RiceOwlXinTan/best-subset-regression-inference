################################################################################
# Linear regression post-selection inference utilities
#
# This script collects the core functions used in the paper experiments and the
# application example. It is designed to be sourced from other scripts, e.g.:
#   source('source code.R')
#
# Main entry points (p-values / intervals):
#   - selective_P_classic(), selective_P(): Monte Carlo selective p-values
#   - selective_P_fast(), selective_P_classic_fast(): optimized variants
#   - saturated_p(): saturated-model based selective p-value via truncated normal
#   - naive_p(): naive (non-selective) p-value
#   - saturated_interval(), naive_interval(): confidence intervals
#
# Distribution helpers:
#   - truncated_cdf_log(): numerically stable CDF for truncated normal on unions
#   - truncated_pdf(), truncated_cdf(): supporting helpers
#
# Notes
# - Interval conventions in this code use a list of numeric length-2 vectors
#   (each interval is c(lower, upper)). Intervals may be unordered and may
#   overlap; some helpers coalesce them internally.
# - This file currently contains sanity-check / plotting code below the
#   distribution helpers. Those lines will execute if you source this file.
#   For non-interactive use, you may want to source selectively or comment out
#   the test block.
################################################################################


###########################
## ---------- Helpers: numeric stability and interval utilities ----------

#' Numerically stable log-sum-exp
#'
#' Computes log(sum(exp(v))) in a stable way.
#'
#' @param v Numeric vector.
#' @return A single numeric value.
logsumexp <- function(v) {
  if (!length(v)) return(-Inf)
  m <- max(v)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(v - m)))
}

#' Numerically stable log(exp(a) - exp(b)) for a >= b
#'
#' @param log_a log(a) where a > 0.
#' @param log_b log(b) where b > 0.
#' @return log(a - b).
logdiffexp <- function(log_a, log_b) {
  if (log_b > log_a) stop("logdiffexp: log_b > log_a")
  log_a + log1p(-exp(log_b - log_a))
}

#' Merge overlapping/touching intervals
#'
#' @param intervals List of numeric vectors c(a,b). Intervals may overlap and
#'   may be unordered.
#' @return A list of disjoint, sorted intervals.
coalesce_intervals <- function(intervals) {
  B <- t(vapply(intervals, function(z) sort(z)[1:2], numeric(2)))
  B <- B[B[,2] > B[,1], , drop = FALSE]
  if (!nrow(B)) return(list())
  B <- B[order(B[,1], B[,2]), , drop = FALSE]
  out <- list(); cur <- B[1,]
  for (i in seq_len(nrow(B))[-1]) {
    if (B[i,1] <= cur[2]) {
      cur[2] <- max(cur[2], B[i,2])
    } else {
      out[[length(out)+1]] <- cur
      cur <- B[i,]
    }
  }
  out[[length(out)+1]] <- cur
  out
}

#' Stable log{ Φ(b) - Φ(a) } for Normal(mu, sigma^2)
#'
#' Uses both lower-tail and upper-tail representations and chooses the more
#' numerically stable one.
#'
#' @param a Lower truncation bound.
#' @param b Upper truncation bound.
#' @param mu Mean.
#' @param sigma Standard deviation (positive).
#' @return log(P(a < Z <= b)) where Z ~ N(mu, sigma^2).
logPhiDiff <- function(a, b, mu, sigma) {
  # assume a < b
  lb <- pnorm(b, mu, sigma, lower.tail = TRUE,  log.p = TRUE)
  la <- pnorm(a, mu, sigma, lower.tail = TRUE,  log.p = TRUE)
  via_lower <- logdiffexp(lb, la)                      # log{Φ(b)-Φ(a)}
  
  ub <- pnorm(b, mu, sigma, lower.tail = FALSE, log.p = TRUE)
  ua <- pnorm(a, mu, sigma, lower.tail = FALSE, log.p = TRUE)
  via_upper <- logdiffexp(ua, ub)                      # log{Φ̄(a)-Φ̄(b)}
  
  # they are equal in exact arithmetic; pick the numerically better one
  if (!is.finite(via_lower)) return(via_upper)
  if (!is.finite(via_upper)) return(via_lower)
  if (via_lower > via_upper) via_lower else via_upper   # closer to 0
}

## ---------- Truncated Normal helpers (log-stable) ----------

#' Truncated normal CDF on a union of intervals (log-stable)
#'
#' Computes the CDF of a Normal(mu, sigma^2) distribution truncated to a union
#' of intervals. The union is specified as a list of c(a,b) bounds and may be
#' unordered and overlapping.
#'
#' @param x Numeric vector of evaluation points.
#' @param mu Mean.
#' @param sigma Standard deviation (positive).
#' @param intervals List of numeric vectors c(lower, upper).
#' @param log.p If TRUE, return log(CDF) rather than CDF.
#' @return Numeric vector of same length as x.


truncated_cdf_log <- function(x, mu, sigma, intervals, log.p = FALSE) {
  if (!is.finite(sigma) || sigma <= 0) stop("sigma must be positive")
  
  intervals <- coalesce_intervals(intervals)
  if (!length(intervals)) {
    warning("No valid intervals.")
    return(if (log.p) rep(-Inf, length(x)) else rep(0, length(x)))
  }
  
  B <- vapply(intervals, identity, numeric(2))
  a <- B[1,]; b <- B[2,]     # sorted, disjoint
  n <- length(a)
  
  # precompute log masses and their prefix log-sum-exp
  log_mass  <- mapply(logPhiDiff, a, b, MoreArgs = list(mu = mu, sigma = sigma))
  log_total <- logsumexp(log_mass)
  if (!is.finite(log_total)) {
    warning("Total truncated mass ≈ 0 under (mu, sigma).")
    return(if (log.p) rep(-Inf, length(x)) else rep(0, length(x)))
  }
  
  prefix_lse <- numeric(n)
  acc <- -Inf
  for (i in seq_len(n)) { acc <- logsumexp(c(acc, log_mass[i])); prefix_lse[i] <- acc }
  
  a_first <- a[1]; b_last <- b[n]
  
  f_one <- function(xi) {
    if (xi <= a_first) return(if (log.p) -Inf else 0)
    if (xi >= b_last)  return(if (log.p) 0    else 1)
    
    # full intervals completely to the left of xi
    k <- findInterval(xi, b, rightmost.closed = TRUE)  # count of b <= xi
    log_left <- if (k > 0) prefix_lse[k] else -Inf
    
    # partial mass if xi is inside interval j
    j <- findInterval(xi, a)                           # largest j with a[j] <= xi
    log_part <- if (j > 0 && xi < b[j]) {
      logPhiDiff(a[j], xi, mu, sigma)
    } else -Inf
    
    log_num <- if (is.finite(log_left) || is.finite(log_part)) {
      logsumexp(c(log_left, log_part))
    } else -Inf
    
    out_log <- log_num - log_total
    if (log.p) out_log else exp(out_log)
  }
  
  vapply(x, f_one, numeric(1))
}


#' Truncated normal CDF (linear-space implementation)
#'
#' This version works in probability space and can lose precision in extreme
#' tails. Prefer truncated_cdf_log() for numerically challenging settings.
#'
#' @inheritParams truncated_cdf_log
#' @return Numeric vector of CDF values.
truncated_cdf <- function(x, mu, sigma, intervals) {
  # Precompute the total probability mass within all intervals
  interval_probs <- sapply(intervals, function(interval) {
    a <- interval[1]
    b <- interval[2]
    pnorm(b, mu, sigma) - pnorm(a, mu, sigma)
  })
  total_mass <- sum(interval_probs)
  
  # Check if total_mass is effectively zero
  if (abs(total_mass) < .Machine$double.eps) {
    warning("The specified intervals contain (almost) no probability mass for the given mean and sigma.")
    # Option 1: Return NaN
    return(rep(0, length(x)))
    # Option 2: Stop execution
    # stop("The specified intervals contain no probability mass.")
  }
  
  # Precompute pnorm values for the interval bounds
  interval_pnorm <- sapply(intervals, function(interval) {
    c(pnorm(interval[1], mu, sigma), pnorm(interval[2], mu, sigma))
  })
  
  # Vectorized CDF computation
  compute_cdf <- function(x_single) {
    if (x_single < intervals[[1]][1]) {
      return(0)
    }
    if (x_single > intervals[[length(intervals)]][2]) {
      return(1)
    }
    # Determine contribution from intervals
    cdf_value <- 0
    for (i in seq_along(intervals)) {
      a <- intervals[[i]][1]
      b <- intervals[[i]][2]
      p_a <- interval_pnorm[1, i]
      p_b <- interval_pnorm[2, i]
      
      if (x_single >= b) {
        # Full interval contributes
        cdf_value <- cdf_value + (p_b - p_a)
      } else if (x_single >= a) {
        # Partial contribution
        cdf_value <- cdf_value + (pnorm(x_single, mu, sigma) - p_a)
        break
      }
    }
    return(cdf_value / total_mass)
  }
  
  # Apply the optimized function to each element of x
  cdf_values <- sapply(x, compute_cdf)
  return(cdf_values)
}



#' Truncated normal PDF on a union of intervals
#'
#' @inheritParams truncated_cdf_log
#' @return Numeric vector of PDF values.
truncated_pdf <- function(x, mu, sigma, intervals) {
  # total mass in allowed intervals
  interval_probs <- sapply(intervals, function(interval) {
    a <- interval[1]; b <- interval[2]
    pnorm(b, mu, sigma) - pnorm(a, mu, sigma)
  })
  total_mass <- sum(interval_probs)
  
  if (abs(total_mass) < .Machine$double.eps) {
    warning("The specified intervals contain (almost) no probability mass for the given mean and sigma.")
    return(rep(0, length(x)))
  }
  
  # indicator: is x inside ANY interval?
  in_any <- sapply(x, function(xx) {
    any(vapply(intervals, function(iv) xx >= iv[1] && xx <= iv[2], logical(1)))
  })
  
  # base normal pdf
  f0 <- dnorm(x, mean = mu, sd = sigma)
  
  # truncated pdf
  out <- ifelse(in_any, f0 / total_mass, 0)
  return(out)
}


#' Reference CDF for disjoint intervals (linear-space)
#'
#' Helper used in the sanity checks.
#'
#' @inheritParams truncated_cdf_log
#' @return Numeric vector of CDF values.
trunc_cdf_linear <- function(x, mu, sigma, intervals) {
  B <- t(vapply(intervals, function(z) sort(z)[1:2], numeric(2)))
  B <- B[order(B[,1]), , drop=FALSE]
  a <- B[,1]; b <- B[,2]
  Fa <- pnorm(a, mu, sigma); Fb <- pnorm(b, mu, sigma)
  mass <- Fb - Fa; total <- sum(mass)
  sapply(x, function(xi) {
    if (xi <= a[1]) return(0)
    if (xi >= b[length(b)]) return(1)
    s <- sum(mass[b <= xi])
    j <- which(a < xi & xi < b)
    if (length(j)) s <- s + (pnorm(xi, mu, sigma) - Fa[j[1]])
    s/total
  })
}


## ---------- Sanity checks (executed on source) ----------
# The following block provides quick numerical checks and plots for the
# truncated normal helpers.

###########test 

set.seed(1)

# 1) Single interval: closed form check
mu <- 0; sigma <- 1
x  <- seq(-10, 10, length.out = 1000)
iv <- list(c(-1, 0.5))
cdf <- truncated_cdf(x, mu, sigma, iv)
ref <- (pnorm(x, mu, sigma) - pnorm(-1, mu, sigma)) /
  (pnorm( 0.5, mu, sigma) - pnorm(-1, mu, sigma))
max(abs(cdf - ref))  # ~ 1e-12 to 1e-15


plot(cdf, type = 'l')
plot(ref, type = 'l')
plot(cdf - ref)


iv2 <- list(c(8, 9),
            c(-9, -8))


cdf2 <- truncated_cdf(x, mu, sigma, iv2)
ref2 <- trunc_cdf_linear(x, mu, sigma, iv2)
max(abs(cdf2 - ref2))  # ~ 1e-12 (linear may lose precision in extreme tails)

plot(x, cdf2, 
     type="l", 
     main="Disjoint intervals [-2,-1] ∪ [1,2]",
     xlab="x", ylab="CDF"); lines(x, ref2, lty=2)


# 3) Tiny tail mass: should be stable
iv3 <- list(c(6, 7))
x3  <- seq(5.5, 7.5, length.out=501)
cdf3 <- truncated_cdf(x3, mu, sigma, iv3)
range(cdf3)  # should go from ~0 to 1 smoothly; no underflow warnings

# 4) Overlapping intervals: coalescing makes result equal to merged interval
iv_overlap <- list(c(-1, 1.5), c(1, 3))    # merge -> [-1, 3]
cdf_overlap <- truncated_cdf(x, mu, sigma, iv_overlap)
cdf_merged  <- truncated_cdf(x, mu, sigma, list(c(-1, 3)))
max(abs(cdf_overlap - cdf_merged))  # 0


plot(x, cdf, type="l", main="Truncated N(0,1) on [-1,2]",
     xlab="x", ylab="CDF"); lines(x, ref, lty=2)


plot(x3, cdf3, type="l", main="Tiny tail mass [6,7]", xlab="x", ylab="CDF")



#' Test whether values fall inside any interval
#'
#' @param x Numeric vector.
#' @param intervals List of numeric vectors c(lower, upper).
#' @return Logical vector of same length as x.
in_any_interval <- function(x, intervals) {
  sapply(x, function(xi) {
    any(vapply(intervals,
               function(iv) iv[1] <= xi && xi <= iv[2],
               logical(1)))
  })
}

# Usage examples:





#' Selective p-value via Monte Carlo (classic implementation)
#'
#' Generates draws from the null (by modifying the target contrast in the
#' estimated coefficient vector) and re-runs best-subset selection each draw.
#' A draw contributes to the selective distribution if it re-selects the
#' originally selected model.
#'
#' This is the most direct version and can be slow for large p or many draws.
#'
#' @param X_dt Design matrix (n x p), without intercept.
#' @param y Response vector (length n) or n x 1 matrix.
#' @param selected Logical vector of length p indicating the selected model.
#' @param alls Logical matrix (S x p); each row is a competing model.
#' @param contrast1 Integer index of the coefficient being tested (among the
#'   selected variables), or 'new' for prediction at a new point.
#' @param contrast2 Optional second index for contrasts of two coefficients.
#' @param new_xpoint Numeric vector of length p (used only when contrast1 == 'new').
#' @param statistic Model selection criterion: 'aic', 'bic', or 'aicc'.
#' @param draws Number of Monte Carlo draws.
#' @param null_value Null value for the tested coefficient/contrast.
#' @param true_sigma Noise level used for simulation.
#' @return Two-sided selective p-value.
selective_P_classic = function(
    X_dt, 
    y, 
    selected, 
    alls, 
    contrast1 = 1,
    contrast2 = NA,
    new_xpoint, 
    statistic='aic', 
    draws = 1000,
    null_value = 0,
    true_sigma = 1){
  
  #input:
  #- X_dt: design matrix
  #- y: response
  #- selected: the selected model (vector: boolean element indicating whether the predictor is included or not)
  #- alls: all other models to be compard with the selected model (binary matrix: each row represent a model)
  #- new_xpoint: new x point
  #- statistic: model selection criteria, must be one of 'aic', 'bic', 'aicc'
  #- alpha: significance level
  #- sigma2: ground truth noise level in the model
  
  #output:
  #- the post-selection confidence interval for predicted mean at new_xpoint
  
  set.seed(1)
  candidates = c()
  everything = c()
  
  
  n = dim(X)[1]
  
  if (contrast1 == 'new'){
    eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  }
  else if( is.na(contrast2) ){
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))

  } 
  else {
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    eta = eta[contrast1+1, ] - eta[contrast2+1, ]
  }
  
  
    cee=eta[contrast1+1,]/c(t(eta[contrast1+1,])%*%eta[contrast1+1,])
    P_select=diag(n)-cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    
    obs = eta[contrast1+1,]%*%y
    
    Eb = eta%*%y
    Eb[contrast1+1,1] = null_value
    
  while  (length(candidates) < draws){
    

    sigma_true = true_sigma
    
    y  = cbind(1,X_dt[,selected])%*%Eb + c(rmvnorm(1, mean = rep(0,n), sigma = diag(n)*sigma_true))
    
    z = y-cee%*%t(eta[contrast1+1,])%*%y
    
    eta_y_sample = eta[contrast1+1,]%*%y
    
    
    df$Y=c(y)
    
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
    
    sample_selected = reg_summary$which[size,-1]
    
    if ( mean(sample_selected == selected) ==1 ){
      candidates = c(candidates, eta_y_sample)
    }
    everything = c(everything,eta_y_sample )
  }
  
  p = 2*min(mean(candidates>obs[1,1]),
            mean(candidates<obs[1,1]))
  
  return (p)
}




#' Selective p-value via interval truncation (classic implementation)
#'
#' Computes a selection region in terms of an excluded set of intervals for the
#' target statistic, then performs Monte Carlo sampling restricted to the
#' allowed region.
#'
#' @inheritParams selective_P_classic
#' @return Two-sided selective p-value.
selective_P=function(
    X_dt, 
    y, 
    selected, 
    alls, 
    contrast1 = 1,
    contrast2 = NA,
    new_xpoint, 
    statistic='aic', 
    draws = 1000,
    null_value = 0,
    true_sigma = 1
    ){
  
  set.seed(1)
  #input:
  #- X_dt: design matrix
  #- y: response
  #- selected: the selected model (vector: boolean element indicating whether the predictor is included or not)
  #- alls: all other models to be compard with the selected model (binary matrix: each row represent a model)
  #- new_xpoint: new x point
  #- statistic: model selection criteria, must be one of 'aic', 'bic', 'aicc'
  #- alpha: significance level
  #- sigma2: ground truth noise level in the model
  
  #output:
  #- the post-selection confidence interval for predicted mean at new_xpoint
  
  
  n=length(y)
  if(statistic=='aic'){
    cons=2
  }else if(statistic=='bic'){
    cons=log(n)
  }else if(statistic=='aicc'){
    # pass
  }else{
    print("must be either aic, aicc or bic!")
    return(NA)
  }
  
  candidates = c()
  everything = c()
  
  if (contrast1 == 'new'){
    eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  }
  else if( is.na(contrast2) ){
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    #eta = eta[contrast1+1, ]
  } 
  else {
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    #eta = eta[contrast1+1, ] - eta[contrast2+1, ]
  }
  
  cee=eta[contrast1+1,]/c(t(eta[contrast1+1,])%*%eta[contrast1+1,])
  P_select=diag(n)-cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
  
  
  obs = eta[contrast1+1,]%*%y
  
  Eb = eta%*%y
  Eb[contrast1+1,1] = null_value
  
  while  (length(candidates) < draws){
    
    #record the excluded intervals for those TruncNormal with 2 intervals
    exclude_lb=rep(NA,nrow(alls))
    exclude_ub=exclude_lb
    
    #record the included intervals for those TruncNormal with 1 interval
    lb_max =-Inf
    ub_min = Inf
    
    sigma_true = true_sigma
    
    y_sample  = cbind(1,X_dt[,selected])%*%Eb + c(rmvnorm(1, mean = rep(0,n), sigma = diag(n)*sigma_true))
    
    # z = rnorm(50, 0, sqrt((1-t(eta)%*%eta))*sigma_true )
    # z = c(rmvnorm(1, mean = rep(0,n), sigma = (diag(n)-cee %*% t(eta))*sigma_true))
    
    z = y_sample - cee%*%t(eta[contrast1+1,])%*%y_sample
    
    eta_y_sample = eta[contrast1+1,]%*%y_sample
    
    for (s in 1:nrow(alls)) {
      #skip comparing with the selected model itself or its superset
      if(all(alls[s,]-selected>=0)){
        next
      }
      
      #exp weight involving cons
      s_hat = sum(selected)#size of selected model
      s_tp  = sum(alls[s,])#size of model alls[s,]
      if(statistic=='aic' | statistic=='bic'){
        wt=exp(cons*(s_hat-s_tp)/n)
      }
      if(statistic=='aicc'){
        wt=exp((2*(s_hat-s_tp)+2*s_hat*(s_hat+1)/(n-s_hat-1)-2*s_tp*(s_tp+1)/(n-s_tp-1))/n)
      }
      #P matrix for model s
      P_s = diag(n)-cbind(1,X_dt[,alls[s,]])%*%solve(t(cbind(1,X_dt[,alls[s,]]))%*%cbind(1,X_dt[,alls[s,]]))%*%t(cbind(1,X_dt[,alls[s,]]))
      #find the coef in quadratic form: A*x^2+B*x+C>0
      Pdiff = P_s-P_select*wt
      A  = c(t(cee)%*%Pdiff%*%cee)
      B1 = c(t(z)%*%Pdiff%*%cee*2)
      C  = c(t(z)%*%Pdiff%*%z)
      #solve the quadratic inequality
      Delta = B1^2 - 4 * A * C
      if(Delta>0){
        r1=(-B1 - sqrt(Delta)) / (2 * A)
        r2=(-B1 + sqrt(Delta)) / (2 * A)
        if(A>0){
          #then there must be 2 intervals as domains, and they have one of {-Inf, Inf} as the tip end
          exclude_lb[s]=min(r1,r2)
          exclude_ub[s]=max(r1,r2)
        }else{
          lb_max=max(min(r1,r2),lb_max)
          ub_min=min(max(r1,r2),ub_min)
        }
      }else{
        if(A<0){
          lb_max=0
          ub_min=0
          return(NA)
        }#otherwise, the interval is the whole R line.
      }
    }
    
    #remove NAs so that what's left always come in pairs
    #the domain is the intersection of {(lb[2i-1], ub[2i-1])union(lb[2i], ub[2i])}
    
    exclude_lb=exclude_lb[!is.na(exclude_lb)]
    exclude_ub=exclude_ub[!is.na(exclude_ub)]
    
    #excluded intervals
    
    if (length(exclude_lb) > 0){
      exclude_intervals=interval_union(Intervals(cbind(exclude_lb,exclude_ub)))
      #also note intervals outside (lb_max, ub_min) are also excluded
      
      if(lb_max>-Inf){
        exclude_intervals=interval_union(exclude_intervals,Intervals(c(-Inf,lb_max)))
      }
      if(ub_min<Inf){
        exclude_intervals=interval_union(exclude_intervals,Intervals(c(ub_min, Inf)))
      }
    }
    else{
      exclude_intervals = Intervals(c(-10000000, -9999999), closed = c(TRUE, TRUE))
    }

    intervals = conv_int(intervals::interval_complement(exclude_intervals) )
    
    if (in_any_interval( eta_y_sample, intervals)){
      candidates = c(candidates, eta_y_sample)
    }
    everything = c(everything, eta_y_sample)
  }
  
  p = 2*min(mean(candidates>obs[1,1]),
            mean(candidates<obs[1,1]))
  return (p)
}


#' Selective p-value via interval truncation (optimized)
#'
#' Optimized version of selective_P() that precomputes QR decompositions for
#' candidate models and uses quadratic-form identities to avoid repeatedly
#' forming projection matrices.
#'
#' @inheritParams selective_P
#' @param CDF If TRUE, return one-sided tail probability P(T >= obs) instead of
#'   the two-sided p-value.
#' @return Two-sided selective p-value (or one-sided tail probability if CDF=TRUE).
selective_P_fast <- function(
    X_dt, 
    y, 
    selected, 
    alls, 
    contrast1 = 1,
    contrast2 = NA,
    new_xpoint, 
    statistic = 'aic', 
    draws = 500,
    null_value = 0,
    true_sigma = 1,
    CDF = FALSE
) {
  set.seed(1)
  
  n <- length(y)
  if (!statistic %in% c("aic","bic","aicc")) stop("statistic must be one of 'aic','bic','aicc'")
  cons <- switch(statistic, aic = 2, bic = log(n), aicc = NA_real_)  # aicc handled per-candidate
  
  if (!requireNamespace("intervals", quietly = TRUE)) {
    stop("Package 'intervals' is required (install.packages('intervals'))")
  }
  
  ## --- Selected model matrices and contrast ---
  Xs <- cbind(1, X_dt[, selected, drop = FALSE])   # n x q
  qsel <- ncol(Xs)
  if (qsel < 2L) stop("selected model must include at least one regressor beyond intercept")
  
  # Build eta exactly as in the slow function, then pick the row (contrast1+1)
  XtXiXt <- solve(t(Xs) %*% Xs) %*% t(Xs)          # q x n
  eta_all <- XtXiXt                                 # rows are contrasts for coefficients (incl. intercept)
  if (identical(contrast1, 'new')) {
    eta_all <- (solve(t(Xs) %*% Xs) %*% t(Xs))     # same as XtXiXt
    eta_row <- (c(1, new_xpoint[selected]) %*% solve(t(Xs) %*% Xs)) %*% t(rep(1, qsel)) # not used; keep API
    # For consistency with slow code:
    eta_all <- (diag(qsel) %*% XtXiXt)
  }
  # pick the contrast row (contrast1+1); slow code uses that row later
  eta_row <- eta_all[contrast1 + 1L, , drop = TRUE]   # length n
  # cee = eta / (eta'eta)
  eta2   <- sum(eta_row * eta_row)
  cee    <- eta_row / eta2
  cee2   <- sum(cee * cee)
  
  # Projection pieces for selected model via QR
  qr_sel    <- qr(Xs)
  Q_sel     <- qr.Q(qr_sel)                       # n x qsel (orthonormal cols)
  qcee_sel  <- as.vector(crossprod(Q_sel, cee))
  qcee_sel2 <- sum(qcee_sel * qcee_sel)
  
  # Useful observed quantities
  obs <- drop(crossprod(eta_row, y))              # eta_row' y
  
  ## --- Candidate models: precompute Q and fixed weights parts ---
  S <- nrow(alls)
  # keep models that are NOT supersets of selected (same skip rule as slow code)
  keep_idx <- which(apply(alls, 1L, function(row) any((row - selected) < 0)))
  
  Q_list    <- vector("list", length(keep_idx))
  qcee_list <- vector("list", length(keep_idx))
  wt        <- numeric(length(keep_idx))
  
  s_hat    <- sum(selected)
  s_sizes  <- rowSums(alls)
  
  for (k in seq_along(keep_idx)) {
    s <- keep_idx[k]
    Xk  <- cbind(1, X_dt[, alls[s,] == 1L, drop = FALSE])
    qrk <- qr(Xk)
    Qk  <- qr.Q(qrk)
    Q_list[[k]]    <- Qk
    qcee_list[[k]] <- as.vector(crossprod(Qk, cee))
    
    if (statistic %in% c("aic","bic")) {
      wt[k] <- exp(cons * (s_hat - s_sizes[s]) / n)
    } else { # aicc
      sh <- s_hat; st <- s_sizes[s]
      wt[k] <- exp((2*(sh - st) + 2*sh*(sh+1)/(n - sh - 1) - 2*st*(st+1)/(n - st - 1)) / n)
    }
  }
  
  ## --- Build Eb exactly like the slow code (replace target with null) ---
  Eb <- as.vector(eta_all %*% y)
  Eb[contrast1 + 1L] <- null_value
  
  ## --- Sampling loop (mirrors slow logic) ---
  candidates <- numeric(0L)
  everything <- numeric(0L)
  
  while (length(candidates) < draws) {
    # y_sample = Xs %*% Eb + epsilon,   epsilon ~ N(0, true_sigma * I)
    y_sample <- drop(Xs %*% Eb) + rnorm(n, mean = 0, sd = sqrt(true_sigma))
    
    # z = y_sample - cee %*% (eta_row' y)
    z <- y_sample - cee * drop(crossprod(eta_row, y_sample))
    
    eta_y_sample <- drop(crossprod(eta_row, y_sample))
    
    # Track excluded intervals (two-interval case) and the single include interval bounds
    excl_lb <- numeric(length(keep_idx))
    excl_ub <- numeric(length(keep_idx))
    fill    <- 0L
    
    lb_max <- -Inf
    ub_min <-  Inf
    
    # Precompute z projections onto selected Q
    qz_sel  <- as.vector(crossprod(Q_sel, z))
    qz_sel2 <- sum(qz_sel * qz_sel)
    z2      <- sum(z * z)
    zcee    <- sum(z * cee)
    
    # Iterate candidate models (fast quadratic via Q)
    for (k in seq_along(keep_idx)) {
      Qk   <- Q_list[[k]]
      qz   <- as.vector(crossprod(Qk, z))
      qz2  <- sum(qz * qz)
      qce  <- qcee_list[[k]]
      qce2 <- sum(qce * qce)
      w    <- wt[k]
      
      # A = cee' (P_s - w P_sel) cee, etc., using P = I - QQ'
      # => A = (1-w)*||cee||^2 - (||Qk'cee||^2 - w||Qsel'cee||^2)
      A  <- (1 - w) * cee2 - (qce2 - w * qcee_sel2)
      # B1 = 2 z' (P_s - w P_sel) cee
      #    = 2[(z'cee - (Qk'z)·(Qk'cee)) - w(z'cee - (Qsel'z)·(Qsel'cee))]
      B1 <- 2 * ( zcee - sum(qz * qce) - w * (zcee - sum(qz_sel * qcee_sel)) )
      # C = z' (P_s - w P_sel) z
      #   = (1-w)||z||^2 - (||Qk'z||^2 - w||Qsel'z||^2)
      C  <- (1 - w) * z2   - (qz2 - w * qz_sel2)
      
      Delta <- B1*B1 - 4*A*C
      if (Delta > 0) {
        r1 <- (-B1 - sqrt(Delta)) / (2*A)
        r2 <- (-B1 + sqrt(Delta)) / (2*A)
        if (A > 0) {
          fill <- fill + 1L
          excl_lb[fill] <- min(r1, r2)
          excl_ub[fill] <- max(r1, r2)
        } else {
          lb_max <- max(lb_max, min(r1, r2))
          ub_min <- min(ub_min, max(r1, r2))
        }
      } else if (A < 0) {
        # empty domain – match slow code behavior
        lb_max <- 0; ub_min <- 0
        return(NA_real_)
      }
    } # candidates
    
    # Build excluded intervals exactly like the slow code
    if (fill > 0L) {
      excl <- intervals::Intervals(cbind(excl_lb[seq_len(fill)], excl_ub[seq_len(fill)]))
      if (is.finite(lb_max)) {
        excl <- intervals::interval_union(excl, intervals::Intervals(c(-Inf, lb_max)))
      }
      if (is.finite(ub_min)) {
        excl <- intervals::interval_union(excl, intervals::Intervals(c(ub_min,  Inf)))
      }
    } else {
      # sentinel tiny exclusion + outside (lb_max, ub_min) if finite
      if (is.finite(lb_max) || is.finite(ub_min)) {
        excl <- intervals::Intervals(c(-Inf, lb_max))
        excl <- intervals::interval_union(excl, intervals::Intervals(c(ub_min, Inf)))
      } else {
        excl <- intervals::Intervals(c(-1e7, -9.999999e6), closed = c(TRUE, TRUE))
      }
    }
    
    intervals_allowed <- conv_int(intervals::interval_complement(excl))
    
    if (in_any_interval(eta_y_sample, intervals_allowed)) {
      candidates <- c(candidates, eta_y_sample)
    }
    everything <- c(everything, eta_y_sample)
  }
  
  p <- 2 * min(mean(candidates > obs), mean(candidates < obs))
  if (CDF){
    return(mean(candidates > obs))
  }
  return(p)
}






#' Sample a random unit vector
#'
#' @param n Dimension.
#' @return Numeric vector of length n with Euclidean norm 1.
sample_unit_vector <- function(n) {
  v <- rnorm(n)          # draw from N(0,1)
  v / sqrt(sum(v^2))     # normalize to unit length
}

# Example:
set.seed(1)
uv = sample_unit_vector(5)*2

t(uv)%*%uv

#' Fit OLS under a single coefficient equality constraint
#'
#' Solves least squares with the constraint beta[j] == c0 by reparameterizing
#' the problem and fitting the remaining coefficients.
#'
#' @param X Design matrix (n x p).
#' @param y Response vector (length n).
#' @param j Which coefficient is fixed (1-based index).
#' @param c0 Fixed value for beta[j].
#' @return List with components: beta (length p) and mu (fitted mean, length n).
fitted_under_beta_j_eq_c <- function(X, y, j, c0) {
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  n <- nrow(X)
  p <- ncol(X)
  
  if (j < 1 || j > p) {
    stop("j must be between 1 and ncol(X)")
  }
  
  # Column j and the remaining columns
  xj <- X[, j, drop = FALSE]
  idx_rest <- setdiff(seq_len(p), j)
  
  # Special case: only one column and it's constrained -> no free coefficients
  if (length(idx_rest) == 0) {
    beta <- rep(NA_real_, p)
    beta[j] <- c0
    mu   <- as.vector(c0 * xj)
    return(list(beta = beta, mu = mu))
  }
  
  X_rest <- X[, idx_rest, drop = FALSE]
  
  # Adjust response for the known contribution c0 * xj
  y_star <- y - c0 * as.numeric(xj)
  
  # OLS of y_star on X_rest
  XtX_rest  <- crossprod(X_rest)         # t(X_rest) %*% X_rest
  Xty_rest  <- crossprod(X_rest, y_star) # t(X_rest) %*% y_star
  beta_rest <- solve(XtX_rest, Xty_rest)
  
  # Full beta under the constraint beta_j = c0
  beta <- numeric(p)
  beta[idx_rest] <- as.numeric(beta_rest)
  beta[j] <- c0
  
  # Fitted mean X %*% beta
  mu <- as.vector(X %*% beta)
  
  list(
    beta = beta,
    mu   = mu
  )
}


#' Selective p-value via Monte Carlo re-selection (fast classic)
#'
#' Faster variant of selective_P_classic() that avoids refitting lm() for each
#' candidate subset size by using RSS values returned by regsubsets().
#'
#' @inheritParams selective_P_classic
#' @param sigmaKnown Logical; whether sigma is treated as known.
#' @param CDF If TRUE, return one-sided tail probability P(T >= obs).
#' @param dist_null If TRUE, return the simulated null draws (diagnostics).
#' @return Two-sided selective p-value (or diagnostics depending on flags).
selective_P_classic_fast <- function(
    X_dt, 
    y, 
    selected, 
    alls, 
    contrast1 = 1,
    actual_p  = 1,
    contrast2 = NA,
    new_xpoint, 
    statistic = 'aic', 
    draws = 1000,
    null_value = 0,
    true_sigma = 1,
    sigmaKnown = TRUE,
    CDF = FALSE,
    dist_null = FALSE
){
  stopifnot(statistic %in% c("aic","bic","aicc"))
  
  set.seed(1)
  
  selected_idx <- which(selected)          # indices of selected vars in the original list
  actual_p <- selected_idx[contrast1]       # original index of the k-th selected

  n <- length(y)
  p <- ncol(X_dt)
  
  # Build data frame & formula ONCE
  df <- as.data.frame(X_dt)
  colnames(df) <- paste0("X", seq_len(p))
  df$Y <- y
  full_fml <- as.formula(paste("Y ~", paste(colnames(df)[1:p], collapse = "+")))
  
  # --- Build eta exactly as in your code ---
  Xs <- cbind(1, X_dt[, selected, drop = FALSE])
  XtX_inv <- solve(t(Xs) %*% Xs)
  XtXiXt  <- XtX_inv %*% t(Xs)   # q x n  (q = 1 + |selected|)
  
  if (identical(contrast1, 'new')) {
    eta <- Xs %*% XtX_inv %*% c(1, new_xpoint[selected])  # n x 1
  } else if (is.na(contrast2)) {
    eta <- XtXiXt                                          # q x n
  } else {
    eta <- XtXiXt
    eta <- eta[contrast1 + 1, , drop = FALSE] - eta[contrast2 + 1, , drop = FALSE]
  }
  
  # Pick the row used for the test (same as your code)
  eta_row <- if (identical(contrast1, 'new')) {
    drop(eta)
  } else if (is.na(contrast2)) {
    eta[contrast1 + 1, ]
  } else {
    drop(eta)
  }
  
  cee  <- eta_row / c(crossprod(eta_row, eta_row))
  obs  <- c(crossprod(eta_row, y))
  
  Eb <- as.vector((if (identical(contrast1, 'new')) 
    (XtX_inv %*% c(1, new_xpoint[selected])) %*% t(Xs)
    else 
      XtXiXt) %*% y)
  Eb[contrast1 + 1] <- null_value
  
  candidates <- numeric(0L)
  everything <- numeric(0L)
  
  # helper to compute AIC/AICc from regsubsets summary
  pick_size_from_ic <- function(reg_sum, n, which_stat) {
    if (which_stat == "bic") {
      return(which.min(reg_sum$bic))
    }
    rss <- reg_sum$rss
    k   <- 1 + rowSums(reg_sum$which[, -1, drop = FALSE])  # +1 intercept
    aic <- n * log(rss / n) + 2 * k
    if (which_stat == "aic")  return(which.min(aic))
    # aicc
    aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
    which.min(aicc)
  }
  
  
  null_y_mean = as.numeric(Xs %*% Eb)
  null_y_mean = fitted_under_beta_j_eq_c(X_dt, 
                                         y, 
                                         j = actual_p, 
                                         c0 = null_value)$mu
  
  
  norm = sqrt(t(y-null_y_mean) %*%(y-null_y_mean))
  
  pre_compute = cee * c(crossprod(eta_row, y))
  
  while (length(candidates) < draws){
    # y_sample = Xs %*% Eb + epsilon, epsilon ~ N(0, true_sigma * I)
    
    if (is.na(true_sigma)){
      y_samp = sample_unit_vector(n)*as.vector(norm) + null_y_mean
      }
    else{
        y_samp <- null_y_mean + rnorm(n, 0, true_sigma)
        }
    
    # keep your definitions
    z <- y_samp - pre_compute
    eta_y_sample <- c(crossprod(eta_row, y_samp))
    
    # Update response only; reuse df/columns/formula
    df$Y <- y_samp
    
    # Best-subset per size via leaps
    models <- leaps::regsubsets(
      full_fml,
      data = df, 
      nvmax = p, 
      method = "exhaustive", 
      nbest = 1,
      intercept = TRUE
    )
    
    reg_summary <- summary(models)
    
    # >>> FAST: no refitting. Use RSS from regsubsets <<<
    size <- pick_size_from_ic(reg_summary, n, statistic)
    
    # Which variables selected in this draw?
    sample_selected <- reg_summary$which[size, -1]
    
    if (mean(sample_selected == selected) == 1) {
      candidates <- c(candidates, eta_y_sample)
    }
    everything <- c(everything, eta_y_sample)
  }
  
  if (dist_null){
    
    return(candidates)
    
  }
  
  pval <- 2 * min(mean(candidates > obs), mean(candidates < obs))
  if (CDF){
    return(mean(candidates > obs))
  }
  pval
}







#' Saturated-model selective p-value via truncated normal
#'
#' Computes the selection region as a union of allowed intervals for the target
#' statistic, then evaluates the truncated normal CDF at the observed value.
#'
#' @inheritParams selective_P_classic
#' @param make_plot If TRUE, produces a diagnostic plot (intended for debugging).
#' @param CDF If TRUE, return CDF value (one-sided) rather than two-sided p-value.
#' @param dist_null If TRUE, return a diagnostic density curve rather than p-value.
#' @return Two-sided selective p-value (or diagnostic outputs depending on flags).
saturated_p=function(
    X_dt, 
    y, 
    selected, 
    alls, 
    contrast1 = 1,
    contrast2 = NA,
    new_xpoint, 
    null_value,
    statistic='aic',
    true_sigma = 1,
    make_plot = FALSE,
    CDF = FALSE,
    dist_null = FALSE
){
  #input:
  #- X_dt: design matrix
  #- y: response
  #- selected: the selected model (vector: boolean element indicating whether the predictor is included or not)
  #- alls: all other models to be compard with the selected model (binary matrix: each row represent a model)
  #- new_xpoint: new x point
  #- statistic: model selection criteria, must be one of 'aic', 'bic', 'aicc'
  #- alpha: significance level
  #- sigma2: ground truth noise level in the model
  
  #output:
  #- the post-selection confidence interval for predicted mean at new_xpoint
  
  n=length(y)
  if(statistic=='aic'){
    cons=2
  }else if(statistic=='bic'){
    cons=log(n)
  }else if(statistic=='aicc'){
    # pass
  }
  else{
    print("must be either aic, aicc or bic!")
    return(NA)
  }
  
  #record the excluded intervals for those TruncNormal with 2 intervals
  exclude_lb=rep(NA,nrow(alls))
  exclude_ub=exclude_lb
  #record the included intervals for those TruncNormal with 1 interval
  lb_max=-Inf
  ub_min=Inf
  
  if (contrast1 == 'new'){
    eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  }
  else if( is.na(contrast2) ){
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    eta = eta[contrast1+1, ]
  } 
  else {
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    eta = eta[contrast1+1, ] - eta[contrast2+1, ]
  }
  
  cee=eta/c(t(eta)%*%eta)
  z=y-cee%*%t(eta)%*%y
  P_select=diag(n)-cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
  
  for (s in 1:nrow(alls)) {
    #skip comparing with the selected model itself or its superset
    if(all(alls[s,]-selected>=0)){
      next
    }
    
    #exp weight involving cons
    s_hat=sum(selected)#size of selected model
    s_tp=sum(alls[s,])#size of model alls[s,]
    if(statistic=='aic' | statistic=='bic'){
      wt=exp(cons*(s_hat-s_tp)/n)
    }
    
    if(statistic=='aicc'){
      wt=exp((2*(s_hat-s_tp)+2*s_hat*(s_hat+1)/(n-s_hat-1)-2*s_tp*(s_tp+1)/(n-s_tp-1))/n)
    }
    
    #P matrix for model s
    P_s=diag(n)-cbind(1,X_dt[,alls[s,]])%*%
      solve(t(cbind(1,X_dt[,alls[s,]]))%*%
              cbind(1,X_dt[,alls[s,]]))%*%
      t(cbind(1,X_dt[,alls[s,]]))
    
    #find the coef in quadratic form: A*x^2+B*x+C>0
    
    Pdiff=P_s-P_select*wt
    A  = c(t(cee)%*%Pdiff%*%cee)
    B1 = c(t(z)%*%Pdiff%*%cee*2)
    C  = c(t(z)%*%Pdiff%*%z)
    #solve the quadratic inequality
    Delta=B1^2 - 4 * A * C
    if(Delta>0){
      r1=(-B1 - sqrt(Delta)) / (2 * A)
      r2=(-B1 + sqrt(Delta)) / (2 * A)
      if(A>0){
        #then there must be 2 intervals as domains, and they have one of {-Inf, Inf} as the tip end
        exclude_lb[s]=min(r1,r2)
        exclude_ub[s]=max(r1,r2)
      }else{
        lb_max=max(min(r1,r2),lb_max)
        ub_min=min(max(r1,r2),ub_min)
      }
    }else{
      if(A<0){
        lb_max=0
        ub_min=0
        return(NA)
      }
      #otherwise, the interval is the whole R line.
    }
  }
  
  #remove NAs so that what's left always come in pairs
  #the domain is the intersection of {(lb[2i-1], ub[2i-1])union(lb[2i], ub[2i])}
  
  exclude_lb=exclude_lb[!is.na(exclude_lb)]
  exclude_ub=exclude_ub[!is.na(exclude_ub)]
  #excluded intervals
  
  # print(Intervals(cbind(exclude_lb,exclude_ub)))

  if (length(exclude_lb) > 0){
    
    exclude_intervals= intervals::interval_union(Intervals(cbind(exclude_lb,exclude_ub)))
    
    #also note intervals outside (lb_max, ub_min) are also excluded
    if(lb_max>-Inf){
      exclude_intervals= intervals::interval_union(exclude_intervals,
                                       Intervals(c(-Inf,lb_max)))
    }
    if(ub_min<Inf){
      exclude_intervals= intervals::interval_union(exclude_intervals,
                                       Intervals(c(ub_min, Inf)))
    }
  }
  else{
    exclude_intervals = Intervals(c(-1000000, -999999), 
                                  closed = c(TRUE, TRUE))
  }
  
  # #if sigmahat is not provided, then we need to estimate
  # #mse-based estiamte from the full model
  # sigmahat1=sigma(lm(y~X_dt))
  # #organic lasso estiamte
  # #sigmahat2=olasso_cv(x = X_dt, y =y)$sig_obj
  # #mse-based estimate from the selected model
  # sigmahat3=sigma(lm(y~X_dt[,selected]))
  # 
  # se1=sqrt(sum(eta^2))*sigmahat1
  # #se2=sqrt(sum(eta^2))*sigmahat2
  # se3=sqrt(sum(eta^2))*sigmahat3
  
  if (is.na(true_sigma)){
    
    s = sigma(lm(y~X_dt[,selected]))
    
  } else if (true_sigma == 'organic'){
    
    s = olasso_cv(x = X_dt, y =y)$sig_obj
    
  } else if (true_sigma == 'full'){
    
    s = sigma(lm(y~X_dt))
    
  } else{
    
    s = true_sigma
    
  }
  
  value=c(t(eta)%*%y)
  se = sqrt(sum(eta^2))*s
  
  #intervals = conv_int((interval_complement(exclude_intervals) - value)/se1 )
  #intervals = conv_int((exclude_intervals - value)/se1 )
  #intervals = conv_int(interval_complement(exclude_intervals) )
  # intervals = conv_int( exclude_intervals)
  # if (intervals[[1]][1] == -10000000){
  #   intervals = list(c(-Inf, Inf))
  # }
  
  df = n - 1 - dim(X_dt)[2]
  #print(intervals)
  #lower_critical <- uniroot(function(x) F_Z(x, 0, 1, intervals, df) - alpha / 2,               lower = -10000, upper = 10000)$root
  #upper_critical <- uniroot(function(x) F_Z(x, 0, 1, intervals, df) - (1 - alpha / 2),         lower = -10000, upper = 10000)$root
  #print("att")
  #print(intervals)
  
  # f1=function(mu){
  #   intervals = conv_int((interval_complement(exclude_intervals)-mu)/se1 )
  #   #print((value-mu)/se1)
  #   #print(mu)
  #   return (F_Z((value-mu)/se1, 0, 1, intervals, df)-alpha/2 )
  # }
  # f2=function(mu){
  #   intervals = conv_int((interval_complement(exclude_intervals)-mu)/se1 )
  #   #print((value-mu)/se1)
  #   #print(F_Z((value-mu)/se1, 0, 1, intervals, df))
  #   #print(mu)
  #   return (F_Z((value-mu)/se1, 0, 1, intervals, df)-1+alpha/2 )
  # }
  # 
  # #print('e1')
  # lb1 = nleqslv(value+qnorm(alpha/2)*se1,   f2)$x
  # ub1 = nleqslv(value+qnorm(1-alpha/2)*se1, f1)$x
  
  intervals = conv_int(intervals::interval_complement(exclude_intervals) )

  # intervals = conv_int((interval_complement(exclude_intervals) - null_value)/se )
  
  if (make_plot){
    x  <- seq(-3, 3, length.out = 000)
    iv <- list(c(-1, 0.5))
    cdf <- truncated_cdf_log(x, null_value, se, intervals)
    plot(cdf, type = 'l', x = x)
    abline(v = value)
    print(intervals)
    print(se)
    print(value)
  }
  
  if (dist_null){
    x  <- seq(-3, 3, length.out = 1000)
    iv <- list(c(-1, 0.5))
    pdf <- truncated_pdf(x, null_value, se, intervals)
    return(pdf)
    
  }
  
  p = truncated_cdf_log(x = value, 
                        null_value,
                        se, 
                        intervals)
  
  if (CDF){
    return (p)
  }
  
  return(2*min(p, 1-p))
}




#' Naive (non-selective) p-value
#'
#' Computes a classical p-value for the target contrast ignoring selection.
#' If true_sigma is NA, uses a t-test with sigma estimated from the selected
#' model; otherwise uses a normal reference.
#'
#' @inheritParams saturated_p
#' @return Two-sided p-value.
naive_p=function(
    X_dt, 
    y, 
    selected, 
    contrast1 = 1,
    contrast2 = NA,
    new_xpoint, 
    null_value,
    statistic='aic',
    true_sigma = 1,
    make_plot = FALSE,
    CDF = FALSE,
    dist_null = FALSE
){
  #input:
  #- X_dt: design matrix
  #- y: response
  #- selected: the selected model (vector: boolean element indicating whether the predictor is included or not)
  #- alls: all other models to be compard with the selected model (binary matrix: each row represent a model)
  #- new_xpoint: new x point
  #- statistic: model selection criteria, must be one of 'aic', 'bic', 'aicc'
  #- alpha: significance level
  #- sigma2: ground truth noise level in the model
  
  #output:
  #- the post-selection confidence interval for predicted mean at new_xpoint
  
  n=length(y)
  
  if (contrast1 == 'new'){
    eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  }
  else if( is.na(contrast2) ){
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    eta = eta[contrast1+1, ]
  } 
  else {
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    eta = eta[contrast1+1, ] - eta[contrast2+1, ]
  }
  
  
  
  if (is.na(true_sigma)){
    s = sigma(lm(y~X_dt[,selected]))
    df = n - 1 - sum(selected)
    
    value=c(t(eta)%*%y)
    se = sqrt(sum(eta^2))*s
    
    # 3) t statistic
    t_stat <- as.numeric((value - null_value) / se)
    
    
    if (dist_null){
      x_l  <- seq(-3, 3, length.out = 1000)
      
      # density of X
      f_x <- function(x, mu, se) {
        (1 / se) * dt((x - mu) / se, df = df)
      }
      
      return(f_x(x_l, null_value, se))
    }
    
    
    
    # 4) p-values
    p_two_sided <- 2 * pt(-abs(t_stat), df)

    
    return(p_two_sided)
    
  } else{
    s = true_sigma
  }
  
  
  value=c(t(eta)%*%y)
  se = sqrt(sum(eta^2))*s
  
  #intervals = conv_int((interval_complement(exclude_intervals) - value)/se1 )
  #intervals = conv_int((exclude_intervals - value)/se1 )
  #intervals = conv_int(interval_complement(exclude_intervals) )
  # intervals = conv_int( exclude_intervals)
  # if (intervals[[1]][1] == -10000000){
  #   intervals = list(c(-Inf, Inf))
  # }
  
  

  p = pnorm(q = value, mean = null_value, sd = se)
  
  
  
  if (CDF){
    return(p)
  }
  
  if (dist_null){
    x_l  <- seq(-3, 3, length.out = 1000)

    # density of X
    f_x <- function(x, mu, se) {
      (1 / se) * dnorm((x - mu) / se,mean = null_value, sd = se)
    }
    
    return(f_x(x_l, null_value, se))
  }
 
  
  
  
  return(2*min(p, 1-p))
  
}




#' Generate an AR(1) correlation matrix
#'
#' @param p Dimension.
#' @param rho AR(1) correlation parameter.
#' @return p x p correlation matrix.
ar1_matrix <- function(p, rho) {
  #input:
  #- p: number of variables
  #- rho: correlation paramter
  
  #output:
  #- a 1st-order auto-regressive correlation matrix
  
  exponent <- abs(matrix(1:p - 1, nrow = p, ncol = p, byrow = TRUE) - (1:p-1))
  rho^exponent
}


#' Convert an intervals::Intervals complement into a list of bounds
#'
#' Many functions in this code expect intervals as a list of c(lower, upper).
#' intervals::interval_complement() returns a different structure; this helper
#' converts it.
#'
#' @param ints An object whose elements are interval bounds as returned by
#'   intervals::interval_complement().
#' @return List of numeric vectors c(lower, upper).
conv_int = function(ints){
  n = length(ints)/2
  result = list()
  for (i in 1:n){
    result[[i]] = c(ints[i][[1]],ints[i][[2]])
    
  }
  return(result)
}





#' Naive confidence interval (ignoring selection)
#'
#' @inheritParams naive_p
#' @param contrast Index of the coefficient/contrast (or 'new' for prediction).
#' @param alpha Significance level.
#' @return Numeric vector c(lower, upper).
naive_interval = function(X_dt,
                          y,
                          selected,
                          contrast,
                          all_compete_models,
                          new_xpoint = new_x,
                          statistic='aic',
                          true_sigma = NA,
                          alpha = 0.05){
  
  
  n=length(y)
  
  if (contrast == 'new'){
    eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  }
  else {
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    eta = eta[contrast+1, ]
  } 
  
  
  
  if (is.na(true_sigma)){
    s = sigma(lm(y~X_dt[,selected]))
    df = n - 1 - sum(selected)
    
    value=c(t(eta)%*%y)
    se = sqrt(sum(eta^2))*s
    
    
  } else{
    s = true_sigma
  }
  
  
  value=c(t(eta)%*%y)
  se = sqrt(sum(eta^2))*s
  
  return(c(value+qnorm(alpha/2)*se, value+qnorm(1-alpha/2)*se))
  
}


#' Saturated-model confidence interval
#'
#' Inverts saturated_p() (via root finding) to produce a two-sided interval.
#'
#' @inheritParams naive_interval
#' @return Numeric vector c(lower, upper).
saturated_interval = function(X_dt,
                              y,
                              selected,
                              contrast,
                              all_compete_models,
                              new_xpoint = new_x,
                              statistic='aic',
                              true_sigma = NA,
                              alpha = 0.05){
  
  
  n=length(y)
  
  if (contrast == 'new'){
    eta=cbind(1,X_dt[,selected])%*%solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%c(1,new_xpoint[selected])
  }
  else {
    eta = solve(t(cbind(1,X_dt[,selected]))%*%cbind(1,X_dt[,selected]))%*%t(cbind(1,X_dt[,selected]))
    eta = eta[contrast+1, ]
  } 
  
  
  
  if (is.na(true_sigma)){
    s = sigma(lm(y~X_dt[,selected]))
    df = n - 1 - sum(selected)
    
    value=c(t(eta)%*%y)
    se = sqrt(sum(eta^2))*s
    
    
  } else{
    s = true_sigma
  }
  
  
  value=c(t(eta)%*%y)
  se = sqrt(sum(eta^2))*s
  
  
  
  
  f1 = function(v){
    saturated_p(X_dt = X_dt,
                y,
                selected = selected,
                contrast1 = contrast,
                contrast2 = NA,
                alls = all_compete_models,
                new_xpoint = new_xpoint,
                statistic='aic',
                null_value = v,
                true_sigma = NA,
                make_plot = FALSE,
                CDF = TRUE)-alpha/2
  }
  
  f2 = function(v){
    saturated_p(X_dt = X_dt,
                y,
                selected = selected,
                contrast1 = contrast,
                contrast2 = NA,
                alls = all_compete_models,
                new_xpoint = new_xpoint,
                statistic='aic',
                null_value = v,
                true_sigma = NA,
                make_plot = FALSE,
                CDF = TRUE)-1+alpha/2
  }
  
  L = nleqslv(value+qnorm(alpha/2)*se,     f2)$x
  U = nleqslv(value+qnorm(1-alpha/2)*se,   f1)$x
  return(c(L,U))
  
}




