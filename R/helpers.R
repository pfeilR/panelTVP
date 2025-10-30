# compute acceptance rates for parameters sampled with Metropolis-Hastings
accept.rate <- function(accept, mcmc.opt){
  accept_burn <- accept[(mcmc.opt$burnin+1):mcmc.opt$chain.length]
  accept_burn_thin <- accept_burn[seq(1,length(accept_burn),by = mcmc.opt$thin)]
  final.accept.rate <- paste("Acceptance Rate =",
                             sum(accept_burn_thin) / length(accept_burn_thin))
  return(final.accept.rate)
}


# preparation of response for R2-step. However, scaling is not done here!
prep_y <- function(y, Time, n){
  yh <- matrix(y, nrow = n, ncol = Time)
  return(yh)
}

# constructZ: construct the regressor matrix for sampling beta and theta in step R2
# Z consists of the stacked matrices Z_t
# beta_tilde is the [Time,d] matrix of noncentered time-varying effect

constructZ <- function(X, Time, timeidx, beta_tilde){

  dimX=dim(X)[2]

  hZ <- list()
  for(it in 1:Time){
    d1 <- X[timeidx == it, ,drop = FALSE]
    if(dimX>1){
      d2 <- diag(beta_tilde[it,], nrow = dimX, ncol = dimX)
    }else{
      d2 <- beta_tilde[it]
    }
    hZ[[it]] <- cbind(d1,d1 %*% d2)
    names(hZ[[it]]) <- NULL
  }
  Z <- do.call("rbind",hZ)
  return(Z)

}

# construct.lp: constructs the linear predictor \eta_{it} =X_{it} beta_t
construct.lp <- function(X, Time, timeidx, betat){
  hlp <- list()
  for(it in 1:Time){
    Xt <- X[timeidx==it, ,drop = FALSE]
    hlp[[it]] <- Xt%*%betat[it,]
    names(hlp[[it]]) <- NULL
  }
  lp <- do.call("rbind",hlp)
  return(lp)
}

#------------------------------------------------------------------------------
# transform_to_centered: from the non-centered to the centered parametrisation
#                        (from beta_tilde to betat)
# beta_tilde is a (T+1)x d matrix
transform_to_centered <- function(beta_tilde, alpha, d){
  TT <- dim(beta_tilde)[1]
  beta <- rep(alpha[1:d], each = TT)
  theta <- alpha[(d+1):(2*d)]
  betat <- matrix(beta, ncol = d) +
    beta_tilde %*% diag(theta, nrow = d, ncol = d)
  return(betat)
}

#--------------------------
# transform_to_noncentered: from the centered to the non-centered parametrisation
#                           (from betat to beta_tilde)
transform_to_noncentered <- function(betat, alpha, d){
  TT <- dim(betat)[1]
  beta <- rep(alpha[1:d], each = TT)
  theta <- alpha[(d+1):(2*d)]
  beta_tilde <- (betat - matrix(beta, ncol = d)) %*%
    diag(1/theta, nrow = d, ncol = d)
  return(beta_tilde)
}

#------------------------------------------------------------------------------
#rGIG_helper: approximates the GIG random variables by gamma and
#inverse gamma for numerically small values of the parameters
#The code is very similar to/translated from the shrinkTVP package
#from Knaus et al. 2019 using version 1.1.1.
rGIG_helper <- function(n, lambda, chi, psi){

  if(chi == 0){
    chi <- .Machine$double.xmin
  }

  if((chi < 11 * .Machine$double.eps)){

    if(lambda > 0){
      res <- rgamma(n, lambda, psi/2)
    }else{
      res <- 1/rgamma(n, -lambda, chi/2)
    }
  }else{
    if(psi < 11 * .Machine$double.eps){
      if(lambda > 0){
        res <- rgamma(n, lambda, psi/2)

      }else{
        res <- 1/rgamma(n, -lambda, chi/2)

      }
    }else{
      res <- GIGrvg::rgig(n = n, lambda = lambda, chi = chi, psi = psi)
    }
  }

  return(res)
}


######
#bessel_k: Function that computes the Bessel Function values
#by either using the base implementation or the Bessel package
#The code is similar to the shrinkTVP package
#from Knaus et al. 2019 using version 1.1.1.

bessel_k <- function(x, nu, bessel_pkg = FALSE){
  if(bessel_pkg){
    res <- Bessel::besselK.nuAsym(x, nu, k.max = 4, expon.scaled = FALSE, log = TRUE)
  }else{
    res <- log(base::besselK(x, nu, expon.scaled = TRUE)) - x
  }
  return(res)
}
