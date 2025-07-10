stepR <- function(response,
                  df,
                  prior.reg,
                  sigma2v = NULL,
                  alpha,
                  W.sparse = NULL,
                  estimation,
                  mcmc.opt,
                  i){

  # Shrinkage prior ------------------------------------------------------------

  if(prior.reg$type %in% c("rw1", "rw2")){

    ## Step DR-1: sample beta_tilde using AWOL ---------------------------------

    if(estimation == "Normal"){
      awol <- AWOL(y = response, X = df$X,
                   Time = df$Tmax, timeidx = df$timeidx, d = df$d,
                   alpha = alpha,
                   sigma2 = sigma2v,
                   hyperpar.c = prior.reg$c,
                   prior = prior.reg$type)
    }
    if(estimation == "PG"){
      awol <- AWOL.PG(z = response, X = df$X, n = df$n,
                      Time = df$Tmax, timeidx = df$timeidx, d = df$d,
                      alpha = alpha,
                      hyperpar.c = prior.reg$c,
                      prior = prior.reg$type,
                      W.sparse = W.sparse)
    }

    if(prior.reg$type == "rw1"){
      beta_tilde <- matrix(awol, nrow = df$Tmax+1, ncol = df$d, byrow = TRUE)
    }

    if(prior.reg$type == "rw2"){
      beta_tilde <- matrix(awol, nrow = df$Tmax, ncol = df$d, byrow = TRUE)
    }

    ## Step DR-2 / TR-2: sample beta and theta ---------------------------------

    A0 <- c(prior.reg$tau, prior.reg$xi)

    if(prior.reg$type == "rw1"){
      Z <- constructZ(df$X, Time = df$Tmax, timeidx = df$timeidx,
                      beta_tilde = beta_tilde[-1,])
    }
    if(prior.reg$type == "rw2"){
      Z <- constructZ(df$X, Time = df$Tmax, timeidx = df$timeidx,
                      beta_tilde = beta_tilde)
    }

    if(estimation == "Normal"){
      alpha <- sample_alpha(y = response, d = df$d, Z = Z, A0 = A0, sigma2 = sigma2v)
    }
    if(estimation == "PG"){
      alpha <- sample_alpha.PG(z = response, d = df$d, Z = Z, A0 = A0, W = W.sparse)
    }

    ## Step DR-3 / TR-3: ancillarity-sufficiency interweaving strategy ---------

    if(mcmc.opt$asis){
      res_ASIS <- ASIS(beta_tilde = beta_tilde,
                       tau = prior.reg$tau,
                       xi = prior.reg$xi,
                       d = df$d,
                       hyp.c = prior.reg$c,
                       alpha = alpha,
                       reg.type = prior.reg$type)
      alpha <- res_ASIS$alpha
      beta_tilde <- res_ASIS$beta_tilde
    }

    ## Step DR-4: sample the hyperparameters -----------------------------------

    ### Step DR-4a: sample a_tau and a_xi (if metropolis == TRUE) --------------

    if(prior.reg$learn.a.xi){

      if(i>1){ # metropolis step starting in second iteration

        ## a_xi

        a_xi_state <- MH_step(a = prior.reg$a.xi,
                              alphapart = alpha[(df$d+1):(2*df$d)],
                              iota = prior.reg$iota.xi,
                              prior_hp1 = prior.reg$nu.xi,
                              prior_hp2 = prior.reg$b.xi*prior.reg$nu.xi,
                              k = prior.reg$kappa.xi,
                              accept = prior.reg$xi.accept,
                              target.rate = prior.reg$target.rate.xi)

        if(a_xi_state[[1]] != prior.reg$a.xi){
          prior.reg$xi.accept[i] <- 1
        } else{
          prior.reg$xi.accept[i] <- 0
        }
        prior.reg$a.xi <- a_xi_state[[1]]
        prior.reg$a.xi[prior.reg$a.xi>10^11]=10^11
        prior.reg$a.xi[prior.reg$a.xi<0.1^7]=0.1^7
        prior.reg$iota.xi <- a_xi_state[[2]]

      }

    }

    if(prior.reg$learn.a.tau){

      if(i>1){ # metropolis step starting in second iteration

        ## a_tau

        a_tau_state <- MH_step(a = prior.reg$a.tau,
                               alphapart = alpha[1:df$d],
                               iota = prior.reg$iota.tau,
                               prior_hp1 = prior.reg$nu.tau,
                               prior_hp2 = prior.reg$b.tau*prior.reg$nu.tau,
                               k = prior.reg$kappa.tau,
                               accept = prior.reg$tau.accept,
                               target.rate = prior.reg$target.rate.tau)

        if(a_tau_state[[1]] != prior.reg$a.tau){
          prior.reg$tau.accept[i] <- 1
        } else{
          prior.reg$tau.accept[i] <- 0
        }
        prior.reg$a.tau <- a_tau_state[[1]]
        prior.reg$a.tau[prior.reg$a.tau>10^11]=10^11
        prior.reg$a.tau[prior.reg$a.tau<0.1^7]=0.1^7
        prior.reg$iota.tau <- a_tau_state[[2]]

      }

    }

    ### Step DR-4b: sample xi, tau ---------------------------------------------

    prior.reg$xi <- sample_GIG(a = prior.reg$a.xi, l = prior.reg$kappa.xi,
                               par = alpha[(df$d+1):(2*df$d)])
    prior.reg$xi[prior.reg$xi>10^11]=10^11
    prior.reg$xi[prior.reg$xi<0.1^15]=0.1^15

    prior.reg$tau <- sample_GIG(a = prior.reg$a.tau, l = prior.reg$kappa.tau,
                                par = alpha[1:df$d])
    prior.reg$tau[prior.reg$tau>10^11]=10^11
    prior.reg$tau[prior.reg$tau<0.1^15]=0.1^15

    ### Step DR-4c: Sample kappa_xi, kappa_tau ---------------------------------

    if(prior.reg$learn.kappa.xi){
      prior.reg$kappa.xi <- sample_G(a = prior.reg$a.xi,
                                     d = df$d,
                                     par = prior.reg$xi,
                                     prior_hp1 = prior.reg$d.xi,
                                     prior_hp2 = prior.reg$e.xi)
      prior.reg$kappa.xi[prior.reg$kappa.xi>10^11]=10^11
      prior.reg$kappa.xi[prior.reg$kappa.xi<0.1^15]=0.1^15
    }

    if(prior.reg$learn.kappa.tau){
      prior.reg$kappa.tau <- sample_G(a = prior.reg$a.tau,
                                      d = df$d,
                                      par = prior.reg$tau,
                                      prior_hp1 = prior.reg$d.tau,
                                      prior_hp2 = prior.reg$e.tau)
      prior.reg$kappa.tau[prior.reg$kappa.tau>10^11]=10^11
      prior.reg$kappa.tau[prior.reg$kappa.tau<0.1^15]=0.1^15
    }

    # transform to centered parameterization for saving results
    betat <- transform_to_centered(beta_tilde = beta_tilde, alpha = alpha, d = df$d)
    if(prior.reg$type == "rw1") betat <- betat[-1,]

    ## Step DR-5 / TR-5: Random sign switch for theta --------------------------

    rsignsw <- base::sample(x = c(-1,1), size = df$d, replace = TRUE)
    theta <- alpha[(df$d+1):(2*df$d)]
    alpha[(df$d+1):(2*df$d)] <- theta*rsignsw

  }

  # Independence prior ---------------------------------------------------------

  else{

    if(estimation == "Normal"){
      betat <- sample.beta.ind(y = response, df = df, sigma2 = sigma2v, B0 = prior.reg$B0)
    }
    if(estimation == "PG"){
      betat <- sample.beta.ind.PG(z = response, df = df, W.sparse = W.sparse, B0 = prior.reg$B0)
    }

  }

  return(list(betat = betat, alpha = alpha, prior.reg = prior.reg))

}

# Additional Functions called by stepR -----------------------------------------

## Step 1 (time-varying effects)

AWOL <- function(y, X, Time, timeidx,d, alpha, sigma2, hyperpar.c, prior){

  #create X_t
  Xt <- lapply(1:Time, FUN = function(t){
    X[timeidx==t, ,drop = FALSE]
  })

  #S_t=(X_t')*(X_t) and is precomputed
  St <- lapply(Xt, FUN = function(x) crossprod(x))

  #diagonal matrix of process standard deviations
  Theta <- diag(alpha[-(1:d)], nrow = d, ncol = d)

  #create Omega using the Matrix package:
  omegaParts <- list()
  Id <- diag(d)

  if(prior == "rw1"){ # Omega_00 exists only when we start at \Tilde{\beta}_0
    #Omega_00:
    omegaParts[[1]] <- (1/hyperpar.c)*Id+Id
    #Omega_tt:
    omegatt <- lapply(1:(Time-1), FUN = function(t){
      omegatt <- (Theta%*%St[[t]]%*%Theta)/sigma2+2*Id
    })
  }

  if(prior == "rw2"){
    #Omega_11:
    omegaParts[[1]] <- (1/hyperpar.c)*Id+Id+(Theta%*%St[[1]]%*%Theta)/sigma2
    #Omega_tt:
    omegatt <- lapply(2:(Time-1), FUN = function(t){
      omegatt <- (Theta%*%St[[t]]%*%Theta)/sigma2+2*Id
    })
  }

  #Omega_TT
  omegaTT <- (Theta%*%St[[Time]]%*%Theta)/sigma2+Id

  omegaParts <- append(omegaParts, omegatt)
  if(prior == "rw1"){
    omegaParts[[Time+1]] <- omegaTT
  }
  if(prior == "rw2"){
    omegaParts[[Time]] <- omegaTT
  }

  # create the tridiagonal block matrix:
  #--- set the block matrices on the diagonal
  om <- Matrix::bdiag(omegaParts)

  #-- set the off-diagonal values: all Omega_t,t+1 are just -I:
  diag(om[,(-1:-d)]) <- -1
  diag(om[(-1:-d),]) <- -1

  #create the c vector:
  h_t <- y - X%*%alpha[1:d]

  c_v <- lapply(1:(Time), FUN = function(t){
    X <- Xt[[t]]
    Ft <- (X%*%Theta)
    (t(Ft))%*%h_t[timeidx==t]/sigma2
  })

  if(prior == "rw1"){
    c_v <- c(rep(0, d), unlist(c_v))
  }
  if(prior == "rw2"){
    c_v <- c(unlist(c_v))
  }

  L <- as.matrix(Matrix::chol(om))

  a <- forwardsolve(t(L), c_v)
  if(prior == "rw1"){
    e <- matrix(rnorm(d*(Time+1)))
  }
  if(prior == "rw2"){
    e <- matrix(rnorm(d*Time))
  }

  #return values:
  backsolve(L, a + e)

}
AWOL.PG <- function(z, X, n, Time, timeidx, d, alpha, W.sparse, hyperpar.c, prior){

  #diagonal matrix of process standard deviations
  Theta <- diag(alpha[-(1:d)], nrow = d, ncol = d)

  #create W matrix -> NEW
  Xt <- lapply(1:(Time), FUN = function(t){
    X[timeidx==t, ,drop = FALSE]
  })

  W.list <- lapply(1:Time, function(t) {
    start_idx <- (t-1)*n+1
    end_idx <- t*n
    return(W.sparse[start_idx:end_idx, start_idx:end_idx, drop = FALSE])
  })

  #create Omega using the Matrix package:
  omegaParts <- list()
  Id <- diag(d)

  if(prior == "rw1"){
    #Omega_00:
    omegaParts[[1]] <- (1/hyperpar.c)*Id+Id
    #Omega_tt -> ADAPTED
    omegatt <- lapply(1:(Time-1), FUN = function(t){
      Wt <- W.list[[t]]
      omegatt <- (Theta %*% t(Xt[[t]]) %*% Wt %*% Xt[[t]] %*% Theta) + 2*Id
    })
  }

  if(prior == "rw2"){
    #Omega_11:
    omegaParts[[1]] <-
      (1/hyperpar.c)*Id+Id+(Theta%*%t(Xt[[1]])%*%W.list[[1]]%*%Xt[[1]]%*%Theta)
    #Omega_tt -> ADAPTED
    omegatt <- lapply(2:(Time-1), FUN = function(t){
      Wt <- W.list[[t]]
      omegatt <- (Theta %*% t(Xt[[t]]) %*% Wt %*% Xt[[t]] %*% Theta) + 2*Id
    })
  }

  #Omega_TT -> ADAPTED
  omegaTT <- (Theta%*%t(Xt[[Time]])%*%W.list[[Time]]%*%Xt[[Time]]%*%Theta)+Id

  omegaParts <- append(omegaParts, omegatt)
  if(prior == "rw1"){
    omegaParts[[Time+1]] <- omegaTT
  }
  if(prior == "rw2"){
    omegaParts[[Time]] <- omegaTT
  }

  # create the tridiagonal block matrix:
  #--- set the block matrices on the diagonal
  om <- Matrix::bdiag(omegaParts)

  #-- set the off-diagonal values: all Omega_t,t+1 are just -I:
  diag(om[,(-1:-d)]) <- -1
  diag(om[(-1:-d),]) <- -1

  #create the c vector -> ADAPTED
  h_t <- z - X %*% alpha[1:d]

  c_v <- lapply(1:Time, FUN = function(t){
    Theta%*%t(Xt[[t]])%*%W.list[[t]]%*%h_t[timeidx==t]
  })

  c_v <- lapply(c_v, function(x) as.matrix(x)) # necessary for converting from sparse to matrix

  if(prior == "rw1"){
    c_v <- c(rep(0, d), unlist(c_v))
  }
  if(prior == "rw2"){
    c_v <- c(unlist(c_v))
  }

  L <- as.matrix(Matrix::chol(om))
  a <- forwardsolve(t(L), c_v)
  if(prior == "rw1"){
    e <- matrix(rnorm(d*(Time+1)))
  }
  if(prior == "rw2"){
    e <- matrix(rnorm(d*Time))
  }

  # added on 16.03.25: remove large objects
  rm(W.list, Xt, omegaParts, om)

  #return values:
  backsolve(L, a + e)

}

## Step 2 (fixed effects)

sample_alpha <- function(y, d, Z, A0, sigma2){

  ZZ <- crossprod(Z,Z)

  # create the square root of A0
  A012 <- diag(sqrt(A0))

  # create the matrix A
  At_star <- (A012%*%ZZ%*%A012)/sigma2 + diag(2*d)

  # chol2inv should be faster
  At <- tryCatch({
    A012 %*% chol2inv(base::chol(At_star)) %*% A012
  }, error = function(e){
    warning("Error in chol2inv: StepR2.R")
    A012 %*% base::solve(At_star) %*% A012
  })

  # create the mean vector a
  a <- (base::tcrossprod(At,Z)%*%y)/sigma2

  mvtnorm::rmvnorm(n=1, mean = a, sigma = At)

}
sample_alpha.PG <- function(z, d, Z, A0, W){

  ZW <- t(Z) %*% W
  ZWZ <- ZW %*% Z

  # create the square root of A0
  A012 <- diag(sqrt(A0))

  ZWZ.A012 <- ZWZ %*% A012
  At_star <- (A012%*%ZWZ.A012) + diag(2*d)

  # chol2inv should be faster
  At <- tryCatch({
    right <- exvatools::multd(chol2inv(base::chol(At_star)), A012)
    A012 %*% right
  }, error = function(e){
    warning("Error in chol2inv: StepR2.R")
    right <- exvatools::multd(base::solve(At_star), A012)
    A012 %*% right
  })

  a <- (At %*% ZW) %*% z
  b <- matrix(MASS::mvrnorm(n = 1, mu = a, Sigma = At), ncol = 1)

  return(b)

}

# Step 3 (ASIS)

ASIS <- function(beta_tilde, d, tau, xi, hyp.c, alpha, reg.type){

  #The actual process stds  (generated in step B for beta)
  theta_old <- alpha[(d+1):(2*d)]

  # R3a
  # interweave into the centered parameterisation and store the signs of theta
  betat <- transform_to_centered(beta_tilde = beta_tilde,
                                 alpha = alpha, d = d)
  TT<- dim(betat)[1]
  theta_signs <- sign(theta_old)

  # R3b:
  #sample beta and theta^2 from the respective conditionals:

  #compute the differences (betat - betat-1) with beta_tilde for numerical stability
  diff1 <- diff(beta_tilde)%*%diag(theta_old, nrow = d, ncol = d)
  diff2 <- beta_tilde[1,]*theta_old

  #resampling:
  theta2 <- sample_theta_squared(diff1 = diff1,
                                 diff2 = diff2,
                                 TT = TT, d = d,
                                 hyp.c = hyp.c,
                                 xi = xi,
                                 reg.type = reg.type)

  beta <- resample_beta(betat = betat,
                        theta2 = theta2,
                        tau = tau,
                        hyp.c = hyp.c,
                        d = d)

  # R3c: generate  the process sds  again and update alpha and beta_tilde
  theta <- theta_signs*sqrt(theta2)
  alpha_new <- c(beta, theta)
  beta_tilde_new <- transform_to_noncentered(betat = betat,
                                             alpha = alpha_new,d = d)

  return(list(alpha = alpha_new,beta_tilde = beta_tilde_new))
}
sample_theta_squared <- function(diff1, diff2, TT, d, hyp.c, xi, reg.type){
  if(reg.type == "rw1") {
    par1 <- -TT/2
  }
  if(reg.type == "rw2"){
    par1 <- -(TT-1)/2
  }
  par2 <- 1/xi
  par3 <- sapply(1:d, FUN = function(x) sum((diff1[,x])^2)+
                   ((diff2[x])^2)/hyp.c)
  sapply(1:d, FUN = function(x) rGIG_helper(n = 1, lambda = par1,
                                            chi = par3[x], psi = par2[x]))
}
resample_beta <- function(betat, theta2, tau, hyp.c, d){

  par1 <- sapply(1:d, function(x) (betat[1,x]*tau[x])/(tau[x]+theta2[x]*hyp.c))
  par2 <- sapply(1:d, function(x) (1/(1/tau[x]+1/(theta2[x]*hyp.c))))

  sapply(1:d, function(x) rnorm(1, par1[x], sd = sqrt(par2[x])))
}

# Step 4 (hyperparameter sampling)

log_acceptance_prob <- function(a_new, a, alphapart, b1, b2, k){

  #HW:  problems for very small alpha
  min.num <-10^(-10)
  if (min(abs(alphapart)<min.num)){
    hsmall <- abs(alphapart)<min.num
    alphapart[hsmall]<-0 # HW
  }
  #HW

  d <- length(alphapart)

  if (sum(alphapart==0)>0){
    alphapart <- alphapart[alphapart != 0]
  }

  logk <- log(k)

  # compute the parameters for the Bessel function:
  par1A <- abs(a_new - 0.5)
  par1B <- abs(a - 0.5)

  par2A <- exp(0.5*log(a_new) + 0.5 * logk + log(abs(alphapart)))
  par2B <- exp(0.5*log(a) + 0.5 * logk + log(abs(alphapart)))

  #--------------------------------------------------
  # Use the approximate or actual Bessel function (as in shrinkTVP)
  bA <- (par2A > 50) | (par1A > 50)
  bB <- (par2B > 50) | (par1B > 50)

  besselA <- sapply(1:d, FUN = function(idx){
    bessel_k(x = par2A[idx], nu = par1A, bessel_pkg = bA[idx])
  }
  )

  besselB <- sapply(1:d, FUN = function(idx){
    bessel_k(x = par2B[idx], nu = par1B, bessel_pkg = bB[idx])
  }
  )

  # compute complete log ratio:
  partA <- (log(a_new) - log(a))*(b1 - 1 + 1+d/4)
  partB <- (a_new - a)*(d*(logk/2) - d*log(2) + sum(log(abs(alphapart))) - b2)
  partC <- (d/2)*(a_new*log(a_new) - a*log(a))
  partD <- -d*(lgamma(a_new + 1) - log(a_new) - lgamma(a+1) + log(a)) # correct, just for numerical stability
  partE <- sum(besselA - besselB)
  res <- partA+partB+partC+partD+partE

  return(res)
}
MH_step <- function(a, alphapart, iota, prior_hp1, prior_hp2, k, accept, target.rate){

  if(!is.na(target.rate)){

    # adaptive step
    accept.rate <- sum(accept) / length(accept)
    adapt.factor <- exp(0.01 * (accept.rate - target.rate))
    iota <- iota * adapt.factor

  } else iota <- iota

  #draw a value from the proposal distribution
  alog <- rnorm(1, mean = log(a), sd = iota)
  proposal <- exp(alog)

  #compute the acceptance prob.:
  acc_prob <- log_acceptance_prob(a_new = proposal,
                                  a = a,
                                  alphapart = alphapart,
                                  b1 = prior_hp1,
                                  b2 = prior_hp2,
                                  k = k)
  u <- runif(1)

  if(log(u) < acc_prob){
    res <- proposal
  }else{
    res <- a
  }
  return(list(res, iota))
}
sample_GIG <- function(a, l, par){

  p1 <- a-0.5
  p2 <- a*l

  par <- par^2
  sapply(par, FUN = function(x)
    rGIG_helper(n = 1, lambda = p1, chi = x, psi = p2))
}
sample_G <- function(a, d, par, prior_hp1, prior_hp2){

  p1 <- prior_hp1+a*d
  p2 <- prior_hp2+(sum(par)*a)/2
  rgamma(1,p1, p2)

}

# independence samplers for betat

sample.beta.ind <- function(y, df, sigma2, B0){

  Tmax <- df$Tmax
  d <- df$d

  B0.inv <- Matrix::Diagonal(Tmax*d, 1/B0)
  X.split <- split(seq_len(nrow(df$X)), df$timeidx)
  X.split <- lapply(X.split, function(rows) df$X[rows, , drop = FALSE])
  X.big <- Matrix::bdiag(X.split)
  Bn <- Matrix::solve(Matrix::crossprod(X.big) / sigma2 + B0.inv)
  bn <- Bn %*% Matrix::t(X.big) %*% y / sigma2

  beta.ind <- MASS::mvrnorm(n = 1, mu = as.numeric(bn), Sigma = as.matrix(Bn))
  beta.ind <- matrix(beta.ind, ncol = d, byrow = TRUE)

  return(beta.ind)

}

sample.beta.ind.PG <- function(z, df, W.sparse, B0){

  Tmax <- df$Tmax
  d <- df$d
  n <- df$n

  B0.inv <- Matrix::Diagonal(Tmax*d,1/B0)
  X.split <- split(seq_len(nrow(df$X)), df$timeidx)
  X.split <- lapply(X.split, function(rows) df$X[rows, , drop = FALSE])
  X.big <- Matrix::bdiag(X.split)
  Bn <- Matrix::solve(Matrix::t(X.big) %*% W.sparse %*% X.big + B0.inv)
  bn <- Bn %*% Matrix::t(X.big) %*% W.sparse %*% z

  beta.ind <- MASS::mvrnorm(n = 1, mu = as.numeric(bn), Sigma = as.matrix(Bn))
  beta.ind <- matrix(beta.ind, ncol = d, byrow = TRUE)

  return(beta.ind)

}
