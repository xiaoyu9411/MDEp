library(meta)
library(estmeansd)
library(pracma)

get.scenario <- function(min.val, q1.val, med.val, q3.val, max.val) {
  if (!(missing(min.val) | missing(q1.val) | missing(med.val) | missing(q3.val)
        | missing(max.val))) {
    if (!any(is.na(c(min.val, q1.val, med.val, q3.val, max.val)))) {
      return("S3")
    }
  }
  if (!(missing(q1.val) | missing(med.val) | missing(q3.val))) {
    if (!any(is.na(c(q1.val, med.val, q3.val)))) {
      return("S2")
    }
  }
  if (!(missing(min.val) | missing(med.val) | missing(max.val))) {
    if (!any(is.na(c(min.val, med.val, max.val)))) {
      return("S1")
    }
  }
  stop("Summary measures not in appropriate form. See documentation for
       appropriate forms.")
}
set.qe.fit.control <- function(quants, n, scenario){
  con <- list()
  if (scenario == "S1" | scenario == "S2") {
    con$norm.mu.bounds <- c(quants[1], quants[3])
    med.val <- quants[2]
    if (min(quants) > 0) {
      con$lnorm.mu.start <- log(med.val)
      con$lnorm.mu.bounds <- c(log(quants[1]), log(quants[3]))
    }
  }
  if (scenario == "S3") {
    con$norm.mu.bounds <- c(quants[2], quants[4])
    med.val <- quants[3]
    if (min(quants) > 0) {
      con$lnorm.mu.start <- log(med.val)
      con$lnorm.mu.bounds <- c(log(quants[2]), log(quants[4]))
    }
  }
  
  con$norm.sigma.bounds <- c(1e-3, 50)
  con$lnorm.sigma.bounds = c(1e-3, 10)
  
  mean.hat <- metaBLUE::Luo.mean(quants, n, scenario)$muhat
  sd.hat <- metaBLUE::Wan.std(quants, n, scenario)$sigmahat
  
  con$norm.mu.start <- mean.hat
  con$norm.sigma.start <- sd.hat
  
  con$lnorm.mu.start <- log(mean.hat / sqrt(1 + (sd.hat/mean.hat)^2))
  con$lnorm.sigma.start <- sqrt(log( 1 + (sd.hat/mean.hat)^2))
  
  con$gamma.shape.start <- mean.hat^2/sd.hat^2
  con$gamma.rate.start <- mean.hat/sd.hat^2
  
  con$beta.shape1.start <- mean.hat *
    (((mean.hat * (1 - mean.hat)) / (sd.hat^2)) - 1)
  con$beta.shape2.start <- con$beta.shape1.start * (1 - mean.hat) / mean.hat
  
  start.val <- (mean.hat/sd.hat)^1.086
  obj.fun <- function(k, mean.hat, sd.hat) {
    temp1 <- gamma((k + 2)/k)
    temp2 <- gamma((k + 1)/k)
    temp3 <- sqrt((temp1 / (temp2^2)) - 1)
    ((sd.hat / mean.hat) - temp3)^2
  }
  no.fit <- function(e) {
    return(list(par = NA))
  }
  get.weibull.start <- tryCatch({
    stats::nlminb(start = start.val, objective = obj.fun, lower = 1e-3,
                  mean.hat = mean.hat, sd.hat = sd.hat)
  },
  error = no.fit,
  warning = no.fit
  )
  if (!is.na(get.weibull.start$par)){
    con$weibull.shape.start <- get.weibull.start$par
  } else {
    con$weibull.shape.start <- start.val
  }
  con$weibull.scale.start <- mean.hat / gamma((con$weibull.shape.start + 1)
                                              / con$weibull.shape.start)
  
  con$gamma.shape.bounds <- c(1e-3, 1e2)
  con$gamma.rate.bounds <- c(1e-3, 1e2)
  con$weibull.shape.bounds <- c(1e-3, 1e2)
  con$weibull.scale.bounds <- c(1e-3, 1e2)
  con$beta.shape1.bounds <- c(10^(-3), 40)
  con$beta.shape2.bounds <- c(10^(-3), 40)
  return(con)
}
get.mean.sd <- function(x, family) {
  if (!(family %in% c("normal", "log-normal", "gamma", "weibull", "beta"))) {
    stop("family must be either normal, log-normal, gamma, Weibull, or beta.")
  }
  if (family == "normal") {
    par <- unname(x$norm.par)
    est.mean <- par[1]
    est.sd <- par[2]
  }
  else if (family == "log-normal") {
    par <- unname(x$lnorm.par)
    est.mean <- exp(par[1] + par[2]^2 / 2)
    est.sd <- sqrt((exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2))
  }
  else if (family == "gamma") {
    par <- unname(x$gamma.par)
    est.mean <- par[1] / par[2]
    est.sd <- sqrt(par[1] / (par[2]^2))
  }
  else if (family == "weibull") {
    par <- unname(x$weibull.par)
    est.mean <- par[2]* gamma(1 + 1/par[1])
    est.sd <- sqrt(par[2]^2 * (gamma(1 + 2 / par[1]) -
                                 (gamma(1 + 1 / par[1]))^2))
  }
  else if (family == "beta") {
    par <- unname(x$beta.par)
    est.mean <- par[1]/(par[1]+par[2])
    est.sd <- sqrt(par[1] * par[2] / ((par[1] + par[2])^2 *
                                        (par[1] + par[2] + 1)))
  }
  return(list(est.mean = est.mean, est.sd = est.sd))
}
get.num.input <- function(min.val, q1.val, med.val, q3.val, max.val, n){
  res <- list()
  if (!missing(min.val)){
    res$min.val <- min.val
  }
  if (!missing(q1.val)){
    res$q1.val <- q1.val
  }
  if (!missing(med.val)){
    res$med.val <- med.val
  }
  if (!missing(q3.val)){
    res$q3.val <- q3.val
  }
  if (!missing(max.val)){
    res$max.val <- max.val
  }
  if (!missing(n)){
    res$n <- n
  }
  return(res)
}

#weighted PE
wpe.fit <- function(min.val, q1.val, med.val, q3.val, max.val, n) {
  
  scenario <- get.scenario(min.val, q1.val, med.val, q3.val, max.val)
  
  if (scenario == "S1") {
    probs <- c(0.625 / n, 0.5, 1 - 0.625 / n)
    quants <- c(min.val, med.val, max.val)
  } else if (scenario == "S2") {
    probs <- c(0.25, 0.5, 0.75)
    quants <- c(q1.val, med.val, q3.val)
  } else {
    probs <- c(0.625 / n, 0.25, 0.5, 0.75, 1 - 0.625 / n)
    quants <- c(min.val, q1.val, med.val, q3.val, max.val)
  }
  
  quants[quants == 0] <- 1e-2
  con <- set.qe.fit.control(quants, n, scenario)
  
  no.fit <- function(e) list(par = NA, value = Inf)
  
  p <- probs
  q <- quants
  
  ## ========= Normal =========
  fit.norm <- tryCatch({
    optim(
      par = c(con$norm.mu.start, con$norm.sigma.start),
      fn = function(par) {
        mu <- par[1]; sigma <- par[2]
        if (sigma <= 0) return(Inf)
        p_pred <- pnorm(q, mu, sigma)
        var <- p_pred * (1 - p_pred)
        sum((p - p_pred)^2 / var)
      },
      method = "L-BFGS-B",
      lower = c(con$norm.mu.bounds[1], con$norm.sigma.bounds[1]),
      upper = c(con$norm.mu.bounds[2], con$norm.sigma.bounds[2])
    )
  }, error = no.fit)
  
  ## ========= Log-normal =========
  fit.lnorm <- tryCatch({
    optim(
      par = c(con$lnorm.mu.start, con$lnorm.sigma.start),
      fn = function(par) {
        mu <- par[1]; sigma <- par[2]
        p_pred <- plnorm(q, mu, sigma)
        var <- p_pred * (1 - p_pred)
        sum((p - p_pred)^2 / var)
      },
      method = "L-BFGS-B",
      lower = c(con$lnorm.mu.bounds[1], con$lnorm.sigma.bounds[1]),
      upper = c(con$lnorm.mu.bounds[2], con$lnorm.sigma.bounds[2])
    )
  }, error = no.fit)
  
  
  ## ========= Weibull =========
  fit.weibull <- tryCatch({
    optim(
      par = c(con$weibull.shape.start, con$weibull.scale.start),
      fn = function(par) {
        shape <- par[1]; scale <- par[2]
        p_pred <- pweibull(q, shape, scale)
        var <- p_pred * (1 - p_pred)
        sum((p - p_pred)^2 / var)
      },
      method = "L-BFGS-B",
      lower = c(con$weibull.shape.bounds[1], con$weibull.scale.bounds[1]),
      upper = c(con$weibull.shape.bounds[2], con$weibull.scale.bounds[2])
    )
  }, error = no.fit)
  
  ## ========= Beta =========
  fit.beta <- tryCatch({
    optim(
      par = c(con$beta.shape1.start, con$beta.shape2.start),
      fn = function(par) {
        shape1 <- par[1]; shape2 <- par[2]
        p_pred <- pbeta(q, shape1, shape2)
        var <- p_pred * (1 - p_pred)
        sum((p - p_pred)^2 / var)
      },
      method = "L-BFGS-B",
      lower = c(con$beta.shape1.bounds[1], con$beta.shape2.bounds[1]),
      upper = c(con$beta.shape1.bounds[2], con$beta.shape2.bounds[2])
    )
  }, error = no.fit)
  
  ## ========= Gamma =========
  fit.gamma <- tryCatch({
    optim(
      par = c(con$gamma.shape.start, con$gamma.rate.start),
      fn = function(par) {
        shape <- par[1]; rate <- par[2]
        p_pred <- pgamma(q, shape, rate)
        var <- p_pred * (1 - p_pred)
        sum((p - p_pred)^2 / var)
      },
      method = "L-BFGS-B",
      lower = c(con$gamma.shape.bounds[1], con$gamma.rate.bounds[1]),
      upper = c(con$gamma.shape.bounds[2], con$gamma.rate.bounds[2])
    )
  }, error = no.fit)
  
  valuesls <- c(
    normal = fit.norm$value,
    `log-normal` = fit.lnorm$value,
    gamma = fit.gamma$value,
    weibull = fit.weibull$value,
    beta = fit.beta$value
  )
  
  output <- list(
    norm.par = fit.norm$par,
    lnorm.par = fit.lnorm$par,
    gamma.par = fit.gamma$par,
    weibull.par = fit.weibull$par,
    beta.par = fit.beta$par,
    valuesls = valuesls,
    scenario = scenario,
    num.input = get.num.input(min.val, q1.val, med.val, q3.val, max.val, n)
  )
  
  class(output) <- "wpe.fit"
  output
}
wpe.mean.sd <- function(min.val, q1.val, med.val, q3.val, max.val, n) {
  
  args <- as.list(environment())
  x <- wpe.fit(min.val, q1.val, med.val, q3.val, max.val, n)
  
  selected.dist <- names(which.min(x$values))
  ests <- get.mean.sd(x, selected.dist)
  
  output <- list(
    est.mean = ests$est.mean,
    est.sd = ests$est.sd,
    selected.dist = selected.dist,
    values = x$values,
    args = args,
    scenario = x$scenario,
    fitted.dists = x
  )
  
  for (dist.name in names(x$values)) {
    ests.all <- get.mean.sd(x, dist.name)
    output[[paste0(dist.name, ".est.mean")]] <- ests.all$est.mean
    output[[paste0(dist.name, ".est.sd")]] <- ests.all$est.sd
  }
  
  class(output) <- "wpe.mean.sd"
  output
}

#MDEp
MDEp.fit<-function(min.val, q1.val, med.val, q3.val, max.val, n)
{
  scenario <- get.scenario(min.val, q1.val, med.val, q3.val, max.val)
  
  if (scenario == "S1") {
    probs <- c(0.625 / n, 0.5, 1 - 0.625 / n)
    quants <- c(min.val, med.val, max.val)
  } else if (scenario == "S2") {
    probs <- c(0.25, 0.5, 0.75)
    quants <- c(q1.val, med.val, q3.val)
  } else if (scenario == "S3") {
    probs <- c(0.625/ n, 0.25, 0.5, 0.75, 1 - 0.625/ n)
    quants <- c(min.val, q1.val, med.val, q3.val, max.val)}
  
  if (min(quants == 0)) {
    quants[quants == 0] <- 10^(-2)
  }
  
  con <-  set.qe.fit.control(quants, n, scenario)
  no.fit <- function(e) {
    return(list(par = NA, value = NA))
  }
  reg <- 1e-6
  fit.norm <- tryCatch({
    N<-n
    start_params <- c(mu = con$norm.mu.start,sigma = con$norm.sigma.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      mu <- params[1]
      sigma <- params[2]
      
      # Ensure sigma is positive
      if (sigma <= 0) return(Inf)
      
      # Calculate predicted probabilities
      p_pred <- pnorm(q, mean = mu, sd = sigma)
      
      
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights 
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i,j]<-p_pred[i]*(1-p_pred[j])
          cov[j,i]<-cov[i,j]
        }
      }
      diag(cov) <- diag(cov) + reg
      
      
      weights <- solve(cov) 
      # Calculate weighted sum of squared residuals 
      residuals <- p - p_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$norm.mu.bounds[1], con$norm.sigma.bounds[1]),
      upper = c(con$norm.mu.bounds[2], con$norm.sigma.bounds[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.n<-qnorm(probs, mean = as.numeric(fit.norm$par[1]), sd = as.numeric(fit.norm$par[2]))
  residual.n<-(quants -pre.n)^2
  value.n<-sum(residual.n)
  
  fit.lnorm <- tryCatch({
    N<-n
    start_params <- c(mu = con$lnorm.mu.start,sigma = con$lnorm.sigma.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      mu <- params[1]
      sigma <- params[2]
      
      # Calculate predicted probabilities
      p_pred <- plnorm(q, meanlog = mu, sdlog = sigma)
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i,j]<-p_pred[i]*(1-p_pred[j])
          cov[j,i]<-cov[i,j]
        }
      }
      diag(cov) <- diag(cov) + reg
      weights <- solve(cov) 
      # Calculate weighted sum of squared residuals 
      residuals <- p - p_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$lnorm.mu.bounds[1], con$lnorm.sigma.bounds[1]),
      upper = c(con$lnorm.mu.bounds[2], con$lnorm.sigma.bounds[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.l<-qlnorm(probs, meanlog = as.numeric(fit.lnorm$par[1]), sdlog = as.numeric(fit.lnorm$par[2]))
  residual.l<-(quants -pre.l)^2
  value.l<-sum(residual.l)
  
  
  fit.weibull <- tryCatch({
    N<-n
    start_params <- c(bsa = con$weibull.shape.start, bsb = con$weibull.scale.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]
      
      # Calculate predicted probabilities
      p_pred <- pweibull(q, shape = bsa, scale = bsb)
      
      
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights 
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i,j]<-p_pred[i]*(1-p_pred[j])
          cov[j,i]<-cov[i,j]
        }
      }
      diag(cov) <- diag(cov) + reg
      weights <- solve(cov) 
      # Calculate weighted sum of squared residuals 
      residuals <- p - p_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$weibull.shape.bounds[1],con$weibull.scale.bounds[1]),
      upper = c(con$weibull.shape.bounds[2],con$weibull.scale.bounds[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.w<-qweibull(probs,shape = as.numeric(fit.weibull$par[1]), scale = as.numeric(fit.weibull$par[2]))
  residual.w<-(quants -pre.w)^2
  value.w<-sum(residual.w)
  
  fit.beta <- tryCatch({
    N<-n
    start_params <- c(con$beta.shape1.start, con$beta.shape2.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]
      
      # Calculate predicted probabilities
      p_pred <- pbeta(q, shape1 = bsa, shape2 = bsb)
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights 
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i,j]<-p_pred[i]*(1-p_pred[j])
          cov[j,i]<-cov[i,j]
        }
      }
      diag(cov) <- diag(cov) + reg
      weights <- solve(cov)  
      # Calculate weighted sum of squared residuals 
      residuals <- p - p_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$beta.shape1.bounds[1],con$beta.shape2.bounds[1]),
      upper = c(con$beta.shape1.bounds[2],con$beta.shape2.bounds[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.b<-qbeta(probs,shape1 = as.numeric(fit.beta$par[1]), shape2 = as.numeric(fit.beta$par[2]))
  residual.b<-(quants -pre.b)^2
  value.b<-sum(residual.b)
  
  
  fit.gamma <- tryCatch({
    N<-n
    start_params <- c(con$gamma.shape.start, con$gamma.rate.start)
    p<-probs
    q<-quants 
    
    objective_function_w <- function(params, q, p, N) {
      bsa <- params[1]
      bsb <- params[2]
      
      # Calculate predicted probabilities
      p_pred <- pgamma(q, shape = bsa, rate = bsb)
      
      
      cov<-matrix(NA,nrow=length(p),ncol=length(p))
      # Calculate variance and weights
      for (i in 1:length(p)){
        for (j in i:length(p)){
          cov[i,j]<-p_pred[i]*(1-p_pred[j])
          cov[j,i]<-cov[i,j]
        }
      }
      diag(cov) <- diag(cov) + reg
      weights <- solve(cov)  
      # Calculate weighted sum of squared residuals
      residuals <- p - p_pred
      weighted_residuals <- t(residuals) %*% weights %*% residuals
      weighted_residuals
    }
    
    result_w <- optim(
      par = start_params,
      fn = objective_function_w,
      q = q,
      p = p,
      N = N,
      method = "L-BFGS-B",
      lower = c(con$gamma.shape.bounds[1], con$gamma.rate.bounds[1]),
      upper = c(con$gamma.shape.bounds[2], con$gamma.rate.bounds[2]),
      hessian = FALSE)
    result_w
  }
  ,   error = no.fit)
  
  pre.g<-qgamma(probs,shape = as.numeric(fit.gamma$par[1]), rate = as.numeric(fit.gamma$par[2]))
  residual.g<-(quants -pre.g)^2
  value.g<-sum(residual.g)
  
  values <- c(value.n, value.l, value.g,
              value.w, value.b)
  
  names(values) <- c("normal", "log-normal", "gamma", "weibull", "beta")
  norm.par <- fit.norm$par
  lnorm.par <- fit.lnorm$par
  gamma.par <- fit.gamma$par
  weibull.par <- fit.weibull$par
  beta.par <- fit.beta$par
  num.input <- get.num.input(min.val, q1.val, med.val, q3.val, max.val, n)
  output <- list(norm.par = norm.par, lnorm.par = lnorm.par,
                 gamma.par = gamma.par, weibull.par = weibull.par,
                 beta.par = beta.par, values = values,
                 num.input = num.input, scenario = scenario)
  class(output) <- "MDE.fit"
  return(output)
  
}
MDEp.mean.sd <- function(min.val, q1.val, med.val, q3.val, max.val, n) {
  args <- as.list(environment())
  x <- MDEp.fit(min.val = min.val, q1.val = q1.val, med.val = med.val,
                q3.val = q3.val, max.val = max.val, n = n)
  selected.dist <- names(which.min(x$values))
  ests <- get.mean.sd(x, selected.dist)
  output <- list(est.mean = ests$est.mean, est.sd = ests$est.sd,
                 selected.dist = selected.dist, values = x$values,args=args,
                 scenario = x$scenario, fitted.dists = x)
  for (dist.name in names(x$values)){
    ests.all <- get.mean.sd(x, dist.name)
    output[paste0(dist.name, '.est.mean')] <- ests.all$est.mean
    output[paste0(dist.name, '.est.sd')] <- ests.all$est.sd
  }
  class(output) <- "MDEp.mean.sd"
  return(output)
}