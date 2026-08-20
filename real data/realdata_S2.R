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
wpe.fit<-function(min.val, q1.val, med.val, q3.val, max.val, n){
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

  # normal
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
          if(i==j){
            cov[i,j]<-p_pred[i]*(1-p_pred[j])
          }
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }}

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

  # log-normal
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
          if(i==j){
            cov[i,j]<-p_pred[i]*(1-p_pred[j])
            cov[j,i]<-cov[i,j]
          }
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }}
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

  # weibull
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
          if(i==j){
            cov[i,j]<-p_pred[i]*(1-p_pred[j])
            cov[j,i]<-cov[i,j]
          }
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }}
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

  # beta
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
          if(i==j){
            cov[i,j]<-p_pred[i]*(1-p_pred[j])
            cov[j,i]<-cov[i,j]
          }
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }}
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

  # gamma
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
          if(i==j){
            cov[i,j]<-p_pred[i]*(1-p_pred[j])
            cov[j,i]<-cov[i,j]
          }
          if(i!=j){
            cov[i,j]<-0
            cov[j,i]<-0
          }
        }}
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
  class(output) <- "wpe.fit"
  return(output)

}
wpe.mean.sd <- function(min.val, q1.val, med.val, q3.val, max.val, n) {
  args <- as.list(environment())
  x <- wpe.fit(min.val = min.val, q1.val = q1.val, med.val = med.val,
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
  class(output) <- "wpe.mean.sd"
  return(output)
}


#MDEp
MDEp.fit<-function(min.val, q1.val, med.val, q3.val, max.val, n){
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

  # normal
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

  # log-normal
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

  # weibull
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

  # beta
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

  # gamma
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



ds<-read.csv("realdata.csv") #Data for PHQ-9 scores

colnames(ds)<-c("Study","a","q1","m","q3","b","n")


dstrue<-read.csv("truevalues.csv")#Include true mean and 95% CI
dstrue$SE<-(dstrue$U-dstrue$L)/(2*1.96)

model<-c()
modelqe<-c()
modelwpe<-c()
modelmdep<-c()

for (i in 1:nrow(ds)){
  dat<-ds[i,2:7]
  n<-as.numeric(dat[6])


  a1<-as.numeric(dat[1])/27
  q11<-as.numeric(dat[2])/27
  m1<-as.numeric(dat[3])/27
  q31<-as.numeric(dat[4])/27
  b1<-as.numeric(dat[5])/27


  ds$luo[i]<- (((0.7+0.39/n)*((q11+q31)/2))+((0.3-0.39/n)*m1))*27

  theta2<-2*qnorm((0.75*n-0.125)/(n+0.25))
  ds$shi[i] <-((q31-q11)/theta2)*27

  ds$mqe[i]<-(as.numeric(qe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[1]))*27
  ds$sdqe[i]<-(as.numeric(qe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[2]))*27
  ss<-summary(qe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n))
  modelqe[i]<-rownames(ss)[1]

  ds$mwpe[i]<-(as.numeric(wpe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[1]))*27
  ds$sdwpe[i]<-(as.numeric(wpe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[2]))*27
  swpe<-wpe.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)
  modelwpe[i]<-swpe$selected.dist


  ds$mmdep[i]<-(as.numeric(MDEp.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[1]))*27
  ds$sdmdep[i]<-(as.numeric(MDEp.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)[2]))*27
  smdep<-MDEp.mean.sd(q1.val = q11,med.val = m1,q3.val = q31,n=n)
  modelmdep[i]<-smdep$selected.dist
}


mymeta1<-metagen(TE=luo,seTE=shi/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta1)

mymeta2<-metagen(TE=mqe,seTE=sdqe/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta2)

mymeta3<-metagen(TE=mwpe,seTE=sdwpe/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta3)

mymeta4<-metagen(TE=mmdep,seTE=sdmdep/sqrt(n),studlab=Study,data=ds,method.tau="REML")
summary(mymeta4)

mymeta5<-metagen(TE=Mean,seTE=SE,studlab=Study,data=dstrue,method.tau="REML")
summary(mymeta5)

#table making
extract_meta_row <- function(meta_obj, method_name) {
  s <- summary(meta_obj)

  data.frame(
    Methods = method_name,
    `Pooled Estimate (95% CI)` = sprintf(
      "%.2f\n(%.2f, %.2f)",
      s$random[1],
      s$lower.random[1],
      s$upper.random[1]
    ),
    tau = sprintf("%.2f", s$tau),
    I2 = sprintf("%.1f%%", s$I2 * 100),
    stringsAsFactors = FALSE
  )
}

table_S2 <- rbind(
  extract_meta_row(mymeta5, "True sample mean/SD"),
  extract_meta_row(mymeta1, "Luo/Shi (Wan)"),
  extract_meta_row(mymeta2, "QE"),
  extract_meta_row(mymeta3, "wPE"),
  extract_meta_row(mymeta4, "MDEp")
)

write.csv(table_S2, "table_S2_2.csv", row.names = FALSE)
