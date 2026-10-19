ll.nb <- function(y, r, psi){

  log_likelihood_i <- lgamma(y + r) - lgamma(r) - lgamma(y + 1) +
    r * log(1 - psi) + y * log(psi)
  log_likelihood <- sum(log_likelihood_i)
  return(log_likelihood)

}

stepRisk <- function(y, miss, eta_nb, eta_logit, r){

  p.risk <- plogis(eta_logit)
  v <- 1 - plogis(eta_nb)

  denom <- 1 - p.risk * (1 - v^r)
  denom[abs(denom) < 1e-12] <- 1e-12

  prob_raw <- (p.risk * v^r) / denom
  prob <- pmin(pmax(prob_raw, 0), 1)
  prob[!is.finite(prob)] <- 0

  risk <- as.logical(ifelse(y == 0 | miss,
                            rbinom(nrow(eta_logit), size = 1, prob = prob),
                            1))

  return(risk)

}

# NEW (added on 19.04.2026): Chinese-Restaurant-Table Distribution Gibbs sampling

sample_china <- function(y, eta, r.old, r.alpha, r.beta){

  n <- length(y)

  # sampling of l_it
  ell <- numeric(n)
  for(i in seq_len(n)){
    if(y[i] > 0){
      prob <- r.old / (r.old + seq_len(y[i]) - 1)
      prob <- pmin(pmax(prob, 0), 1)
      ell[i] <- sum(rbinom(y[i],1,prob))
    }
  }

  # sampling of r
  r <- rgamma(1, r.alpha + sum(ell), r.beta + sum(log1p(exp(eta))))

  return(r)

}

# NEW (added on 12.05.2026): Blocked MH-Sampling of (r, beta1)

sample_r_beta1 <- function(y,
                           eta,
                           r,
                           beta1,
                           eps.sd.blocked,
                           r.alpha,
                           r.beta,
                           tau1,
                           accept,
                           target.rate,
                           adapt){

  if(adapt && !is.na(target.rate)){
    accept.rate <- sum(accept) / length(accept)
    adapt.factor <- exp(0.01 * (accept.rate - target.rate))
    eps.sd.blocked <- eps.sd.blocked * adapt.factor
  }

  # old values
  log.r_old <- log(r)
  r_old <- r
  beta1_old <- beta1
  psi_old <- plogis(eta)
  psi_old <- pmin(pmax(psi_old, 1e-10), 1 - 1e-10)

  # candidates
  eps <- rnorm(1, sd = eps.sd.blocked)
  log.r_star <- log.r_old + eps
  r_star <- exp(log.r_star)
  beta1_star <- beta1 - eps
  psi_star <- plogis(eta - eps)
  psi_star <- pmin(pmax(psi_star, 1e-10), 1 - 1e-10)

  # computing acceptance probability
  num <- ll.nb(y = y, r = r_star, psi = psi_star) +
    dgamma(r_star, shape = r.alpha, rate = r.beta, log = TRUE) +
    dnorm(beta1_star, mean = 0, sd = sqrt(tau1), log = TRUE) + log.r_star

  den <- ll.nb(y = y, r = r_old, psi = psi_old) +
    dgamma(r_old, shape = r.alpha, rate = r.beta, log = TRUE) +
    dnorm(beta1_old, mean = 0, sd = sqrt(tau1), log = TRUE) + log.r_old

  alpha <- min(0, num - den)

  u <- runif(1)

  if(alpha < log(u)){
    return(list(r = r_old, beta1 = beta1_old, eps = 0, eps.sd.blocked = eps.sd.blocked))
  } else{
    return(list(r = r_star, beta1 = beta1_star, eps = eps, eps.sd.blocked = eps.sd.blocked))
  }

}

# NEW (added on 19.05.2026): Instead of CRT and Slice, Metropolis for log(r)

sample_r_metro <- function(y, eta, r.old, eps.sd.r, r.alpha, r.beta, accept, target.rate, adapt){

  # adaptive step
  if(adapt && !is.na(target.rate)){
    accept.rate <- sum(accept) / length(accept)
    adapt.factor <- exp(0.01 * (accept.rate - target.rate))
    eps.sd.r <- eps.sd.r * adapt.factor
  }

  rho.old <- log(r.old)

  eps <- rnorm(1, mean = 0, sd = eps.sd.r)

  rho.star <- rho.old + eps
  r.star <- exp(rho.star)

  psi <- pmin(pmax(plogis(eta), 1e-10), 1 - 1e-10)

  num <- ll.nb(y = y, r = r.star, psi = psi) +
    dgamma(r.star, shape = r.alpha, rate = r.beta, log = TRUE) +
    rho.star

  den <- ll.nb(y = y, r = r.old, psi = psi) +
    dgamma(r.old, shape = r.alpha, rate = r.beta, log = TRUE) +
    rho.old

  logacc <- min(0, num - den)

  if(log(runif(1)) < logacc){
    return(list(r = r.star, eps.sd.r = eps.sd.r))
  } else {
    return(list(r = r.old, eps.sd.r = eps.sd.r))
  }

}
