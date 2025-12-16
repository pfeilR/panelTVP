slice_IV_2D <- function(response, residual.stage1, sigma2, rho, prior.var_stage2, prior.rho){

  # transform to make support of both parameters unbounded
  val <- c(log(sigma2), atanh(rho)) # val[1] = log(sigma2), val[2] = Fisher-Z(rho)

  # auxiliary function to compute log-posterior in two dimensions
  log_post_2D <- function(vals){

    sigma2.val <- exp(vals[1])
    rho.val <- tanh(vals[2])
    if(sigma2.val <= 0 || abs(rho.val) >= 1) return(-Inf)
    e <- response - rho.val * sqrt(sigma2.val) * residual.stage1
    log.lik <- sum(dnorm(e, mean = 0, sd = sqrt(sigma2.val * (1 - rho.val^2)), log = TRUE))
    log.prior_sigma <- dgamma(sigma2.val, prior.var_stage2$alpha.sigma,
                              prior.var_stage2$beta.sigma, log = TRUE)
    log.prior_rho <- dnorm(vals[2], mean = prior.rho$mean.rho, sd = prior.rho$sd.rho, log = TRUE)

    return(log.lik + log.prior_sigma + vals[1] + log.prior_rho) # Jacobian for log(sigma^2)

  }

  log.y0 <- log_post_2D(val)
  slice.level <- log.y0 + log(runif(1))
  w <- c(prior.var_stage2$width, prior.rho$width)
  steps <- c(prior.var_stage2$expansion.steps, prior.rho$expansion.steps)
  L <- val - runif(2, 0, w)
  R <- L + w

  # stepping out phase (in 2D!)
  for(d in 1:2){
    J <- floor(runif(1, 0, steps[d]))
    K <- steps[d] - 1 - J
    while(J > 0 && log_post_2D(replace(val, d, L[d])) > slice.level){
      L[d] <- L[d] - w[d]
      J <- J - 1
    }
    while(K > 0 && log_post_2D(replace(val, d, R[d])) > slice.level){
      R[d] <- R[d] + w[d]
      K <- K - 1
    }
  }

  # shrinkage based on 2D-rectangle area
  repeat{
    val.new <- runif(2, L, R)
    if(log_post_2D(val.new) >= slice.level) break
    for(d in 1:2){
      if(val.new[d] < val[d]) L[d] <- val.new[d] else R[d] <- val.new[d]
    }
  }

  # deliver back-transformed value
  sigma2.new <- exp(val.new[1])
  rho.new <- tanh(val.new[2])

  return(list(sigma2 = sigma2.new, rho = rho.new))

}

