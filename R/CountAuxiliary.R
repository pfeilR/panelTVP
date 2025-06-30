NB.para <- function(y,
                    eta = NULL,
                    iota.r = NULL,
                    r.old = NULL,
                    omega.old = NULL,
                    sample.r = FALSE,
                    sample.omega = FALSE,
                    compute.z = FALSE,
                    r.a = NULL,
                    r.b = NULL,
                    r.accept = NULL,
                    r.target.rate = NULL,
                    slice = FALSE,
                    r.alpha = NULL,
                    r.beta = NULL,
                    expansion.steps = NULL,
                    width = NULL,
                    p.overrelax = NULL,
                    accuracy.overrelax = NULL){

  nn <- length(y)

  if(sample.r){

    psi <- 1 / (1 + exp(-eta)) # DON'T USE exp(eta)/(1+exp(eta)) as this is instable!!!
    # psi should not be exactly 0 or 1
    psi <- ifelse(psi == 0, psi + 0.00001, psi)
    psi <- ifelse(psi == 1, psi - 0.00001, psi)
    if(!slice){ # Metropolis-Hastings update
      r.next <- MH.r(y = y, psi = psi, r = r.old, iota.r = iota.r, r.a = r.a, r.b = r.b,
                     r.target.rate = r.target.rate, r.accept = r.accept)
    } else{ # Slice-Sampler update
      r.next <- slice.r(y = y, psi = psi, r = r.old, a = r.alpha, b = r.beta,
                        steps = expansion.steps, w = width, p.overrelax = p.overrelax,
                        acc = accuracy.overrelax)
    }
    return(r.next)

  }

  if(sample.omega){

    omega <- BayesLogit::rpg(num = nn, h = as.numeric(y + r.old), z = eta)
    return(omega)

  }

  if(compute.z){

    z <- (y - r.old)/(2*omega.old)
    return(z)

  }

}

ll.nb <- function(y, r, psi){

  log_likelihood_i <- lgamma(y + r) - lgamma(r) - lgamma(y + 1) +
    r * log(1 - psi) + y * log(psi)
  log_likelihood <- sum(log_likelihood_i)
  return(log_likelihood)

}

MH.r <- function(y, r, psi, iota.r, r.a, r.b, r.target.rate = NULL, r.accept = NULL){

  # optional: adaptive proposal choice
  if(!is.null(r.accept)){
    r.accept.rate <- sum(r.accept) / length(r.accept) # computing acceptance rate
    adapt.factor <- exp(0.01 * (r.accept.rate - r.target.rate))
    iota.r <- iota.r * adapt.factor
  }

  # optional part ends here ---
  r.star <- truncnorm::rtruncnorm(n = 1, a = r.a, b = r.b, mean = r, sd = iota.r)
  num <- ll.nb(y = y, psi = psi, r = r.star) + log(truncnorm::dtruncnorm(r, a = r.a, b = r.b, mean = r.star, sd = iota.r))
  den <- ll.nb(y = y, psi = psi, r = r) + log(truncnorm::dtruncnorm(r.star, a = r.a, b = r.b, mean = r, sd = iota.r))
  log.alpha <- num - den
  if(is.nan(log.alpha)){ # reject if results are very unusual to prevent error
    return(list(r = r, iota.r = iota.r))
  }
  alpha <- exp(min(0, log.alpha))
  u <- runif(1)
  if(alpha > u){
    return(list(r = r.star, iota.r = iota.r))
  } else{
    return(list(r = r, iota.r = iota.r))
  }

}

slice.r <- function(y, psi, r, a, b, steps, w, p.overrelax, acc){

  # auxiliary function to compute (un-normalized) log-posterior under a G(a,b) prior
  log.post.r <- function(log.r.val){
    r.val <- exp(log.r.val)
    if(r.val <= 0) return(-Inf)
    log.prior <- dgamma(r.val, shape = a, rate = b, log = TRUE)
    log.lik <- ll.nb(y = y, psi = psi, r = r.val)
    val <- log.lik + log.prior
    if(!is.finite(val)) return(-Inf)
    return(val)
  }

  # Slice phase
  log.r <- log(r)
  log.y0 <- log.post.r(log.r)
  slice.level <- log.y0 + log(runif(1)) # this gives me the horizontal slice

  # Stepping out phase
  u <- runif(1, max = w) # interval of size w is randomly around log.r
  L <- log.r - u
  R <- log.r + (w - u)
  J <- floor(runif(1, 0, steps)) # to avoid expanding for far too long
  K <- steps - 1 - J
  while(J > 0 && log.post.r(L) > slice.level){
    L <- L - w
    J <- J - 1
  }
  while(K > 0 && log.post.r(R) > slice.level){
    R <- R + w
    K <- K - 1
  }

  overrelax <- sample(c(FALSE, TRUE), size = 1, prob = c(1-p.overrelax, p.overrelax))

  if(overrelax){

    # Overrelaxation phase

    # Step 1: narrow interval until mid-point is inside the slice (or accuracy limit is hit)
    if((R-L) < 1.1 * w){
      repeat{
        mid <- (L+R)/2
        if(acc == 0 || log.post.r(mid) > slice.level) break
        if(log.r > mid){
          L <- mid
        } else{
          R <- mid
        }
        acc <- acc-1
        w <- w/2
      }
    }

    # Step 2: refine endpoints using bisection
    L.hat <- L
    R.hat <- R
    while(acc > 0){
      acc <- acc - 1
      w <- w/2
      if(slice.level >= log.post.r(L.hat+w)){
        L.hat <- L.hat + w
      }
      if(slice.level >= log.post.r(R.hat-w)){
        R.hat <- R.hat - w
      }
    }

    # Step 3: find a suitable candidate by flipping and possibly rejecting
    log.r.new <- L.hat + R.hat - log.r
    if(L > log.r.new || R < log.r.new || slice.level >= log.post.r(log.r.new)){
      log.r.new <- log.r # rejection (necessary for detailed balance)
    }

  } else{

    # Shrinkage phase
    repeat{
      log.r.new <- runif(1, L, R)
      if(log.post.r(log.r.new) >= slice.level) break
      if(log.r.new < log.r){
        L <- log.r.new
      } else{
        R <- log.r.new
      }
    }

  }

  return(list(r = exp(log.r.new)))

}

stepRisk <- function(y, miss, eta_nb, eta_logit, r){

  p.risk <- plogis(eta_logit)
  v <- 1 - plogis(eta_nb)
  risk <- as.logical(ifelse(y == 0 | miss,
                            rbinom(nrow(eta_logit), size = 1, prob =
                                     pmax(0, pmin(1, (p.risk * v^r) / (1 - p.risk * (1 - v^r))))), 1))

  return(risk)

}
