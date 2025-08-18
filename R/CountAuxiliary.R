NB.para <- function(y,
                    eta = NULL,
                    r.old = NULL,
                    omega.old = NULL,
                    sample.r = FALSE,
                    sample.omega = FALSE,
                    compute.z = FALSE,
                    r.alpha = NULL,
                    r.beta = NULL,
                    expansion.steps = NULL,
                    width = NULL,
                    p.overrelax = NULL,
                    accuracy.overrelax = NULL){

  nn <- length(y)

  if(sample.r){

    psi <- 1 / (1 + exp(-eta))
    psi <- ifelse(psi == 0, psi + 0.00001, psi)
    psi <- ifelse(psi == 1, psi - 0.00001, psi)
    r.next <- slice_r(y = y, psi = psi, r = r.old, a = r.alpha, b = r.beta,
                      steps = expansion.steps, w = width, p.overrelax = p.overrelax,
                      acc = accuracy.overrelax)
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

slice_r <- function(y, psi, r, a, b, steps, w, p.overrelax, acc){

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
        acc <- acc - 1
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

  relevant.obs <- (y == 0 | miss)
  p.risk <- pmax(1e-4, pmin(1-1e-4, plogis(eta_logit[relevant.obs])))
  v <- pmax(1e-4, pmin(1-1e-4, 1 - plogis(eta_nb[relevant.obs])))
  risk <- as.logical(ifelse(relevant.obs,
                            rbinom(length(eta_logit[relevant.obs]),
                                   size = 1,
                                   prob = (p.risk * v^r) / (1 - p.risk * (1 - v^r))
                                   ),
                            1)
                     )

  return(risk)

}
