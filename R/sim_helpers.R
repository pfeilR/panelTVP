sim_Gaussian_Probit_Logit_NegBin <- function(n,
                                             Tmax,
                                             beta,
                                             theta,
                                             lambda,
                                             psi,
                                             model,
                                             r = NULL,
                                             sigma2 = NULL){

  t <- Tmax

  # regression effects ---------------------------------------------------------

  d <- length(beta)
  beta0nc <- rnorm(d)
  betatnc <- apply(matrix(c(beta0nc,rnorm(t*d)),byrow=TRUE,ncol=d),2,cumsum)
  D <- diag(theta)
  beta.th <- matrix(beta, nrow = t+1, ncol = d, byrow = TRUE) +
    t(apply(betatnc, 1, function(x) D%*%x))
  betat <- beta.th[-1,]

  # factor loadings ------------------------------------------------------------

  lambdat <- matrix(lambda + psi * rnorm(t), ncol = 1)

  # data generation ------------------------------------------------------------

  indivID <- rep(1:n, times = t)
  W <- matrix(rnorm(n*t*d, mean = 0), nrow = n*t)
  W <- model.matrix(~W[,-1])
  linpred <- matrix(NA,n*t)
  for(tt in 1:t){
    ind <- ((tt-1)*n+1):(n*tt)
    linpred[ind,] <- W[ind,] %*% matrix(betat[tt,])
  }
  if(model == "Gaussian"){
    Sigma <- diag(rep(sigma2, t))
    error <- matrix(mvtnorm::rmvnorm(n, mean = rep(0, t), sigma = Sigma))
  }
  y <- matrix(NA,n*t)
  fi <- base::scale(rnorm(n))
  Fmat <- matrix(NA,n*t)
  for(tt in (1:t)){
    Fmat[((tt-1)*n+1):(n*tt),] <- fi %*% matrix(lambdat[tt,])
  }

  if(model != "Gaussian"){

    eta <- linpred + Fmat

    if(model == "Logit"){
      prob.binary <- plogis(eta)
      y <- rbinom(length(prob.binary), size = 1, prob = prob.binary)
    }

    else if(model == "Probit"){
      prob.binary <- pnorm(eta)
      y <- rbinom(length(prob.binary), size = 1, prob = prob.binary)
    }

    else if(model == "NegBin"){
      mu <- r * exp(eta)
      y <- MASS::rnegbin(length(mu), mu = mu, theta = r)
    }

  } else{

    y <- linpred + Fmat + error

  }

  # Setting up return object ---------------------------------------------------

  W <- W[,-1]
  namess <- c()
  for(i in 1:(d-1)){
    namess[i] <- paste0("W",i)
  }
  if((d-1)==1) W <- as.matrix(W)
  colnames(W) <- namess
  observed <- data.frame(y = y, W, t = rep(1:t, each = n), id = indivID)

  if(model == "Gaussian"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                sigma2 = sigma2,
                model = "Gaussian")

  }

  else if(model == "Probit"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                model = "Probit")

  }

  else if(model == "Logit"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                model = "Logit")

  }

  else if(model == "NegBin"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                r = r,
                model = "NegBin")

  }

  return(ret)

}


sim_ZINB <- function(n,
                     Tmax,
                     beta.nb,
                     theta.nb,
                     lambda.nb,
                     psi.nb,
                     beta.logit,
                     theta.logit,
                     lambda.logit,
                     psi.logit,
                     r){


  t <- Tmax

  # simulation of betat for logit component ------------------------------------

  d.logit <- length(beta.logit)
  beta0nc.logit <- rnorm(d.logit)
  betatnc.logit <- apply(matrix(c(beta0nc.logit, rnorm(t*d.logit)), byrow = TRUE,
                                ncol = d.logit), 2, cumsum)
  D.logit <- diag(theta.logit)
  beta.th.logit <- matrix(beta.logit, nrow = t+1, ncol = d.logit, byrow = TRUE) +
    t(apply(betatnc.logit, 1, function(x) D.logit%*%x))
  betat.logit <- beta.th.logit[-1,]

  # simulation of betat for negative binomial component ------------------------

  d.nb <- length(beta.nb)
  beta0nc.nb <- rnorm(d.nb)
  betatnc.nb <- apply(matrix(c(beta0nc.nb, rnorm(t*d.nb)), byrow = TRUE,
                             ncol = d.nb), 2, cumsum)
  D.nb <- diag(theta.nb)
  beta.th.nb <- matrix(beta.nb, nrow = t+1, ncol = d.nb, byrow = TRUE) +
    t(apply(betatnc.nb, 1, function(x) D.nb%*%x))
  betat.nb <- beta.th.nb[-1,]

  # factor loadings for logit component ----------------------------------------

  lambdat.logit <- matrix(lambda.logit + psi.logit * rnorm(t), ncol = 1)

  # factor loadings for negative binomial component ----------------------------

  lambdat.nb <- matrix(lambda.nb + psi.nb * rnorm(t), ncol = 1)

  # data generation ------------------------------------------------------------

  indivID <- rep(1:n, times = t)

  W.logit <- matrix(rnorm(n*t*d.logit, mean = 0), nrow = n*t)
  W.logit <- model.matrix(~W.logit[,-1])
  W.nb <- matrix(rnorm(n*t*d.nb, mean = 0), nrow = n*t)
  W.nb <- model.matrix(~W.nb[,-1])
  linpred.logit <- matrix(NA,n*t)
  linpred.nb <- matrix(NA,n*t)
  fi.logit <- base::scale(rnorm(n))
  Fmat.logit <- matrix(NA,n*t)
  fi.nb <- base::scale(rnorm(n))
  Fmat.nb <- matrix(NA,n*t)
  for(tt in 1:t){
    ind <- ((tt-1)*n+1):(n*tt)
    linpred.logit[ind,] <- W.logit[ind,] %*% matrix(betat.logit[tt,])
    linpred.nb[ind,] <- W.nb[ind,] %*% matrix(betat.nb[tt,])
    Fmat.logit[ind,] <- fi.logit %*% matrix(lambdat.logit[tt,])
    Fmat.nb[ind,] <- fi.nb %*% matrix(lambdat.nb[tt,])
  }
  eta.logit <- linpred.logit + Fmat.logit
  eta.nb <- linpred.nb + Fmat.nb

  p_at.risk <- plogis(eta.logit)
  w <- rbinom(n*t, size = 1, prob = p_at.risk)
  mu <- r * exp(eta.nb)
  y <- ifelse(w == 1, MASS::rnegbin(length(mu), mu = mu, theta = r), 0)

  # setting up return object ---------------------------------------------------

  W.logit <- W.logit[,-1]
  names.logit <- c()
  for(i in 1:(d.logit-1)){
    names.logit[i] <- paste0("W",i,".logit")
  }
  if((d.logit-1)==1) W.logit <- as.matrix(W.logit)
  colnames(W.logit) <- names.logit
  W.nb <- W.nb[,-1]
  names.nb <- c()
  for(i in 1:(d.nb-1)){
    names.nb[i] <- paste0("W",i,".nb")
  }
  if((d.nb-1)==1) W.nb <- as.matrix(W.nb)
  colnames(W.nb) <- names.nb
  observed <- data.frame(y = y, W.logit, W.nb, t = rep(1:t, each = n), id = indivID)
  ret <- list(observed = observed,
              beta_zinb.inflation = betat.logit,
              beta_zinb.count = betat.nb,
              lambda_zinb.inflation = lambdat.logit,
              lambda_zinb.count = lambdat.nb,
              r = r,
              model = "ZINB")
  return(ret)

}

