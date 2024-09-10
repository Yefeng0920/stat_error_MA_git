# power (two-tail) for meta-analysis
powerMA <- function(mu, SE, alpha = 0.05) {
  1 - pnorm(qnorm(1 - 0.05 / 2) - abs(mu) / SE) + pnorm(-qnorm(1 - 0.05 / 2)-abs(mu) / SE)
} 

# S error for meta-analysis
errorS <- function(mu, se, alpha = 0.05){
  p.u <- 1 - pnorm(qnorm(1 - alpha/2) - abs(mu)/se) 
  p.l <- pnorm(-qnorm(1 - alpha/2) - abs(mu)/se) 
  power <- p.u + p.l 
  errorS <- p.l/power 
  return(errorS)
} 

errorM <- function(mu, se, alpha = 0.05, N = 1e5) {
  est.random <- rnorm(n=N, mean = mu, sd = se)
  sig.index <- suppressWarnings(abs(est.random) > se*qnorm(1 - alpha/2))
  overestimate <- mean(abs(est.random)[sig.index])/abs(mu) 
  absolute_error <- overestimate*abs(mu) - abs(mu)
  relative_error <- absolute_error/(overestimate*abs(mu))
  return(abs(overestimate) %>% round(3))
}



# meta-analysis of magnitude
## folded effect size
folded_es <-function(mean, variance){ # the sampling variance of magnitude   
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_mu
}
## folded error
folded_error <- function(mean, variance){ # the sampling variance of magnitude   
  mu <- mean
  sigma <- sqrt(variance)
  fold_mu <- sigma*sqrt(2/pi)*exp((-mu^2)/(2*sigma^2)) + mu*(1 - 2*pnorm(-mu/sigma))
  fold_se <- sqrt(mu^2 + sigma^2 - fold_mu^2)
  # adding se to make bigger mean
  fold_v <- fold_se^2
  fold_v
}

# functions for the mixture distribution
dmix = function(x,p,m,s){
  drop(p %*% sapply(x, function(x) dnorm(x,mean=m,sd=s)))
}

rmix = function(n,p,m,s){    # sample from a normal mixture
  d=rmultinom(n,1,p)
  rnorm(n,m%*%d,s%*%d)
}

pmix = function(x,p,m,s){ # cumulative distr (vector x)
  drop(p %*% sapply(x, function(x) pnorm(x,mean=m,sd=s)))
}

# p(SNR | z) when snr ~ dmix(p,m,s)
posterior <- function(z,p,m,s) { 
  s2=s # we estimate s of snr directly
  p=p*dnorm(z,m,s2)  
  p <- p/sum(p)         # conditional mixing probs
  pm <- z*s2/(s2+1) + m/(s2+1) # conditional means
  pv <- s2/(s2+1)          # conditional variances
  ps <- sqrt(pv)            # conditional std devs
  data.frame(p,pm,ps)
}

# function also accepts censoring below c1 and above c2.
# uses base R function `constrOptim` to run constrained optimization. 
# we define the log likelihood function and set up constraints such
# that the mixture proportions are non-negative and add up to 1.
mix = function(z, k = 3, c1 = 0, c2 = 10^6, weights = 1){
  # log likelihood function
  loglik = function(theta, z, k, weights = 1){
    p = c(theta[1:(k-1)], 1 - sum(theta[1:(k-1)]))
    s = theta[k:(2*k-1)]
    m = rep(0, k)
    lik1 = (abs(z) < c1) * (pmix(c1, p, m = m, s = s) - pmix(-c1, p, m = m, s = s))
    lik2 = (abs(z) >= c1) * (abs(z) < c2) * dmix(z, p, m = m, s = s)
    lik3 = (abs(z) >= c2) * (pmix(-c2, p, m = m, s = s) + 1 - pmix(c2, p, m = m, s = s))
    lik = lik1 + lik2 + lik3
    return(-sum(weights * log(lik)))   # *minus* the weighted log lik
  }
  
  # set up constraints for optimization
  # The feasible region is defined by ui %*% par - ci >= 0
  ui = c(rep(-1, (k-1)), rep(0, k))         # (k-1) mixture props sum to < 1
  ui = rbind(ui, diag(2*k-1))
  ci = c(-1, rep(0, k-1), rep(0, k))        # Removed constraint about variances being at least 1
  
  # set starting value
  theta0 = c(rep(1/k, (k-1)), c(1.2, 2:k))
  opt = constrOptim(theta = theta0, f = loglik, ui = ui, ci = ci,
                    method = "Nelder-Mead",
                    z = z, weights = weights, k = k,
                    control = list(maxit = 10^4))
  
  # collect the results
  p = c(opt$par[1:(k-1)], 1 - sum(opt$par[1:(k-1)]))  # mixture proportions
  sigma = opt$par[k:(2*k-1)]                          # mixture sds
  m = rep(0, k)                                        # mixture means
  df = data.frame(p = p, m = m, sigma = sigma)
  return(df)
}

# probability of replication when the distribution is asymmetric (non zero-mean)
replcalc2 <- function(zscore,p,m,sigma){ # compute replication probability (predictive power) when original study produced z
  z=abs(zscore)
  pr=dmix(z,p,m,sigma) / (dmix(z,p,m,sigma) + dmix(-z,p,m,sigma)) # pr(z >0 | |z|)
  pr=drop(pr)
  post=posterior(z,p,m,sigma) # p(SNR|z= |z|)
  pm=post$pm
  ps=post$ps
  powpos=1 - pmix(1.96,p=post$p,m=pm,s=sqrt(ps^2 + 1)) # signif given z=|z|
  post=posterior(-z,p,m,sigma) # p(SNR|z=-|z|)
  pm=post$pm
  ps=post$ps
  powneg=pmix(-1.96,p=post$p,m=pm,s=sqrt(ps^2 + 1))   # signif given z=-|z|
  return(as.numeric(pr*powpos + (1-pr)*powneg)) # signif given |z|
}

# probability of replication when the distribution is symmetric (zero-mean)
replcalc <- function(z,p,m,s,multiplier=1){
  post=posterior(abs(z),p=p,m=m,s=s)
  pp=post$p
  pm=sqrt(multiplier)*post$pm
  ps=sqrt(multiplier)*post$ps
  1 - pmix(1.96,p=pp,m=pm,s=sqrt(ps^2 + 1))
}



# odds ratio and base rate
calcRodds <- function(FDR, pwr, alpha = 0.05) {
  R = alpha * (1 - FDR) / (pwr * FDR)
  Pr = R / (R + 1)
  return(Pr)
}
