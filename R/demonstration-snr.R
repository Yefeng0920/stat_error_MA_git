library(flexmix)
library(DEoptim)

### custom functions for fiting mixture models of absolute z-statistics
# these functions are very similar to z-curve
# the main difference is that these models:
# - use all test statistics (while z-curve uses only significant results by default)
# - estimate location and standard deviations of the mixture (while z-curve fixes means to 0:7 and standard deviations to 1)
zdist_pdf   <- function(x, mean, sd){
  l1 <- dnorm(x,  mean, sd)
  l2 <- dnorm(x, -mean, sd)
  return(l1 + l2)
}
fit_wzdist  <- function(x, w, fix_sd = FALSE, min_sd = FALSE){

  if(isFALSE(fix_sd)){
    fit <- DEoptim(
      fn  = function(par, x, w){
        -sum(w*log(zdist_pdf(x, mean = par[1], sd = par[2])))
      },
      x = x,
      w = w,
      lower   = c(0,   if(isFALSE(min_sd)) 0.001 else min_sd),
      upper   = c(100, 100),
      control = list(
        trace = FALSE
      )
    )
    param <- as.list(fit$optim$bestmem)
    names(param) <- c("mean", "sd")
  }else{
    fit <- DEoptim(
      fn  = function(par, x, w, fix_sd){
        -sum(w*log(zdist_pdf(x, mean = par[1], sd = fix_sd)))
      },
      x      = x,
      fix_sd = fix_sd,
      w      = w,
      lower   = c(0),
      upper   = c(100),
      control = list(
        trace = FALSE
      )
    )
    param <- as.list(fit$optim$bestmem)
    names(param) <- c("mean")
  }

  return(param)
}
zdist_clust <- function(formula = .~., fix_sd = FALSE, min_sd = FALSE){

  retval <- new("FLXMC", weighted = TRUE,
                formula = formula, dist = "zdist",
                name = "folded normal distributions")

  retval@defineComponent <- function(param) {

    logLik <- function(x, y) {
      log(zdist_pdf(y, mean = param$mean, sd = param$sd))
    }

    predict <- function(x) {
      matrix(param$mean, nrow = nrow(x), ncol = length(param$mean), byrow = TRUE)
    }

    new("FLXcomponent",
        parameters = param,
        df         = param$df,
        logLik     = logLik,
        predict    = predict)
  }

  retval@fit <- function(x, y, w, ...) {
    param <- fit_wzdist(y, w, fix_sd, min_sd)
    if(!isFALSE(fix_sd)){
      param$sd <- fix_sd
    }
    retval@defineComponent(c(param, df = if(isFALSE(fix_sd)) 2 else 1))
  }
  return(retval)
}
density_zdist_clust <- function(x, fit){
  weights <- fit@size / sum(fit@size)
  param   <- parameters(fit)
  x_dens  <- do.call(rbind, lapply(seq_along(weights), function(i)
    weights[i] * zdist_pdf(x, param["mean",i], param["sd",i])
  ))
  return(colSums(x_dens))
}
density_snr_clust   <- function(x, fit, ...){

  weights <- fit@size / sum(fit@size)
  param   <- parameters(fit)

  # add continuous density
  x_dens  <- do.call(rbind, lapply(seq_along(weights), function(i)
    if(param["sd",i] > 1){
      weights[i] * zdist_pdf(x, param["mean",i], sqrt(param["sd",i]^2 - 1))
    }else{
      rep(0, length(x))
    }
  ))
  if(any(colSums(x_dens) > 0)){
    lines(x, colSums(x_dens), ...)
  }

  # add fixed density
  for(i in seq_along(weights)){
    if(param["sd",i] == 1){
      abline(v = param["mean",i], ..., lwd = 10*weights[i])
    }
  }
  return(invisible(NA))
}

### simulate srn -> z-statistics
# I simulate data from a bivariate distribution of signal-to-noise ratio
# i.e., some (almost null) studies and some considerably powered studies
set.seed(1)
z_srn <- abs(c(rnorm(4000, 0, 0.4), rnorm(2000, 2, 0.1)))
z     <- abs(rnorm(length(z_srn), z_srn, 1))

# check distribution of observed z-statistics and signal-to-noise ratios
# notice that the z-statistics are decreasing despite the signal-to-noise ratio is increasing
# (van Zwet model assumes that signal-to-noise is decreasing from 0
#  --- because the normal distributions are zero centered)
par(mfrow = c(1, 2))
hist(z,     freq = FALSE, breaks = 50, xlab = "z-statistic")
hist(z_srn, freq = FALSE, breaks = 50, xlab = "signal-to-noise ratio")
dev.off()

### fit mixture models that allow for non-centered signal-to-noise (takes a while)
# the first one fixes the standard deviation to 1 (i.e., a z-curve model with variable means)
# the second one fixes the standard deviation to 1 (i.e., a snr-curve model with variable means)
fit1 <- flexmix(z ~ 1, k = 2, model = zdist_clust(fix_sd = 1))
fit2 <- flexmix(z ~ 1, k = 2, model = zdist_clust(min_sd = 1))

summary(fit1)
parameters(fit1)

summary(fit2)
parameters(fit2)

# see the marginal distribution of observed statistics
# notice that both models recovery it almost exactly
hist(z, freq = FALSE, breaks = 50)
curve(density_zdist_clust(x, fit1), add = TRUE, col = "red")
curve(density_zdist_clust(x, fit2), add = TRUE, col = "blue")

# however, they imply vastly different distribution of signal-to-noise ratios
# the first model assumes that all snr are fixed at the estimated means (sd = 1)
# --- therefore you see only the vertical lines
# while the second model assumes that the snr are normally distributed around the estimated means
# --- therefore you can see the distribution
hist(z_srn, freq = FALSE, breaks = 50)
density_snr_clust(seq(0, 3, 0.01), fit1, col = "red")
density_snr_clust(seq(0, 3, 0.01), fit2, col = "blue")
# importantly, both models are wrong despite having the correct number of components
# (and the second one even the correct form)
# and both of them imply very different conclusions about the distribution of signal-to-noise ratios
# --- and the subsequent metrics

# you can try fitting and visualizing the centered model to the data
# the problem is that it forces the mean of the SNR distribution to be centered at zero
# --- i.e., it is impossible to find that there is more powered than non-power studies
# it is a fundamental flaw of the model which hardcodes this assumption
