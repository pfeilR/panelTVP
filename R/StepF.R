stepF <- function(response,
                  df,
                  sigma2v = NULL,
                  W.sparse = NULL,
                  W.dense = NULL,
                  lambda,
                  alpha_lambda,
                  prior.load,
                  estimation){

  res.mat <- matrix(response, byrow = F, nrow = df$n, ncol = df$Tmax)

  # Step F1: sample the subject specific factors -------------------------------

  if(estimation == "Normal"){
    fi <- sample.facs(yit = res.mat, lambda = lambda, sigma2 = sigma2v)
  }
  if(estimation == "PG"){
    fi <- sample.facs.PG(zit = res.mat, lambda = lambda, W.dense = W.dense)
  }

  # Step F2: sample the factor loadings and hyperparameters --------------------

  if(estimation == "Normal"){
    load.draws <- sample_lambda(y = response, fi = fi, Time = df$Tmax,
                                timeidx = df$timeidx, sigma2 = sigma2v,
                                alpha = alpha_lambda, pri.lambda = prior.load)
  }
  if(estimation == "PG"){
    load.draws <- sample_lambda.PG(z = response, zit = res.mat, fi = fi,
                                   Time = df$Tmax, timeidx = df$timeidx,
                                   alpha = alpha_lambda, pri.lambda = prior.load,
                                   W.sparse = W.sparse, W.dense = W.dense)
    rm(W.sparse, W.dense)
  }

  lambda <- load.draws$lambda
  alpha_lambda <- load.draws$alpha
  prior.load <- load.draws$pri.lambda

  # Step F3: Sign-Switch

  # signsw <- base::sample(x = c(-1,1), size = 1, replace = TRUE)
  # lambda <- lambda*signsw
  # fi <- fi*signsw
  # alpha_lambda <- alpha_lambda*signsw

  # Step F4: Identification ----------------------------------------------------

  if(lambda[1] < 0){
    lambda <- -lambda
    fi <- -fi
    alpha_lambda <- -alpha_lambda
  }

  return(list(fi = fi, lambda = lambda, alpha_lambda = alpha_lambda, prior.load = prior.load))

}

# Additional Functions called by stepF -----------------------------------------

## Step 1 (sample factor scores)

# NB: For the ZINB model the factor scores are sampled using the information of
# the data in the at-risk set of the current iteration. However, if a subject
# has repeated measurements that are only not-at risk, then the corresponding
# factor score is sampled from the N(0,1) prior (as suggested in Neelon, 2019).

sample.facs <- function(yit, lambda, sigma2){

  n <- dim(yit)[1]
  Tmax <- dim(yit)[2]
  if(length(lambda)==1){
    Mni <- 1/(1+Tmax*lambda^2/sigma2)
    yh <- lambda*rowSums(yit)/sigma2
  } else{
    Mni <- 1/(1+sum(lambda^2)/sigma2)
    yh <- yit%*%matrix(lambda/sigma2, ncol=1)
  }
  mni <- Mni*yh
  fi <- rnorm(n, mean = mni, sd = sqrt(Mni))
  fi <- fi-mean(fi)
  return(fi)

}
sample.facs.PG <- function(zit, lambda, W.dense){

  n <- dim(zit)[1]
  Tmax <- dim(zit)[2]
  if(length(lambda)==1){
    Mni <- 1/(1+lambda^2*rowSums(W.dense))
    yh <- lambda*rowSums(zit*W.dense)
  } else{
    Mni <- 1/(1+(W.dense%*%matrix(lambda, ncol=1)^2))
    yh <- (zit*W.dense)%*%matrix(lambda, ncol=1)
  }
  mni <- Mni*yh
  fi <- rnorm(n, mean = mni, sd = sqrt(Mni))
  fi <- fi-mean(fi)
  return(fi)

}

## Step 2 (sample factor loadings and hyperparameters)

sample_lambda <- function (y, fi, Time, timeidx, alpha, pri.lambda, sigma2){

  n <- length(fi)

  if(pri.lambda$type=="cps"|pri.lambda$type=="ind"){

    sum.f2 <- sum(fi^2)

    XSX <- sum.f2*diag(1/sigma2, Time) # here diag necessary for dimension
    yitS <- matrix(y,ncol=Time)/sigma2

    if(pri.lambda$type=="cps"){
      XSX <- sum(XSX)
      inv.L0 <- 1/as.vector(pri.lambda$L0[1])
      Ln <- 1/(XSX+inv.L0)
      ln <- Ln*sum(yitS*fi) # equivalent to yitS + matrix(rep(fi,Time ))
      lambda <- rnorm(1, mean=ln, sd=sqrt(Ln))
    }else{ # independence prior
      Ln <- solve(XSX + diag(1/pri.lambda$L0, ncol(XSX)))
      ln <- Ln%*%colSums(yitS*fi)
      lambda <- MASS::mvrnorm(1, mu=ln, Sigma=Ln)
    }

  }

  if(pri.lambda$type=="rw1"|pri.lambda$type=="rw2"){ # random walk priors

    ystar <- y-rep(alpha[1]*fi, Time)
    psi <- alpha[2]

    lambda_tilde <- AWOL_fac(yf = ystar, fi = fi, Time = Time,
                             psi = psi, pri.type = pri.lambda$type,
                             hyp.c = pri.lambda$c,
                             sgma2 = sigma2)

    lambda_tilde <- matrix(lambda_tilde, ncol = 1)

    # sample lambda and psi2
    A0lambda <- c(pri.lambda$phi, pri.lambda$zeta)

    Xf <- matrix(rep(fi, Time),ncol=1)
    if(pri.lambda$type=="rw1"){
      lth <- lambda_tilde[-1]
    } else{
      lth <- lambda_tilde
    }

    Zlambda <- constructZ(X = Xf, Time = Time, timeidx = timeidx, beta_tilde=lth)

    alpha_lambda <- sample_alpha(y = y, d = 1, Z = Zlambda, A0 = A0lambda,
                                 sigma2 = sigma2)

    lambdah <- transform_to_centered(beta_tilde = lambda_tilde,
                                     alpha = alpha_lambda,
                                     d = 1)
    if(pri.lambda$type=="rw1"){
      lambda <- lambdah[1+(1:Time)] # factor loading in observations
    }else{
      lambda <- lambdah
    }

    # sample phi_lambda and zeta_lambda

    phi_lambda <- sample_GIG(a = pri.lambda$a.phi,
                             l = pri.lambda$kappa.phi,
                             par = alpha_lambda[1])
    phi_lambda[phi_lambda>10^11]=10^11
    phi_lambda[phi_lambda<0.1^15]=0.1^15

    zeta_lambda <- sample_GIG(a = pri.lambda$a.zeta,
                              l = pri.lambda$kappa.zeta,
                              par = alpha_lambda[2])
    zeta_lambda[zeta_lambda>10^11]=10^11
    zeta_lambda[zeta_lambda<0.1^15]=0.1^15

    pri.lambda$phi <- phi_lambda
    pri.lambda$zeta <- zeta_lambda

    # sampling of kappa_phi and kappa_zeta

    if(pri.lambda$learn.kappa.phi){
      pri.lambda$kappa.phi <- sample_G(a = pri.lambda$a.phi,
                                       d = 1,
                                       par = pri.lambda$phi,
                                       prior_hp1 = pri.lambda$d.phi,
                                       prior_hp2 = pri.lambda$e.phi)
    }

    if(pri.lambda$learn.kappa.zeta){
      pri.lambda$kappa.zeta <- sample_G(a = pri.lambda$a.zeta,
                                        d = 1,
                                        par = pri.lambda$zeta,
                                        prior_hp1 = pri.lambda$d.zeta,
                                        prior_hp2 = pri.lambda$e.zeta)
    }

  }

  if(pri.lambda$type == "rw1" | pri.lambda$type == "rw2"){
    res <- list(lambda=lambda, alpha=alpha_lambda, pri.lambda=pri.lambda)
  }
  else{
    res <- list(lambda=lambda, alpha=alpha, pri.lambda=pri.lambda)
  }
  return(res)
}
AWOL_fac <- function(yf, fi, Time, psi, pri.type, hyp.c, sgma2){

  n <- length(fi)
  Sf <- rep(psi^2*sum(fi^2),Time)/sgma2

  # construct prior matrix
  if(pri.type=="rw1"){
    d <- Time+1
    h <- list(rep(-1,d-1),
              c(1+1/hyp.c, rep(2, d-2),1))

    L0.inv <- as.matrix(Matrix::bandSparse(d, k=(-1):0, diag=h,symm=TRUE))
    Omega <- diag(c(0,Sf))+L0.inv
  }

  if(pri.type=="rw2"){
    d <- Time
    h <- list(rep(-1,d-1),
              c(1+1/hyp.c, rep(2, d-2),1))
    L0.inv <- as.matrix(Matrix::bandSparse(d, k=(-1):0, diag=h,symm=TRUE))
    Omega <- diag(Sf)+L0.inv
  }

  #create the c vector
  yh <- prep_y(yf, Time, n)
  cv <- matrix(psi*crossprod(fi,matrix(yh, nrow=n, ncol=Time))/sgma2, ncol=1)
  if(pri.type=="rw1"){cv=c(0,cv)}

  Oinv <- solve(Omega)
  m <- Oinv%*%cv
  lambda_tilde <- mvtnorm::rmvnorm(n=1, mean = m, sigma = Oinv)

  return(lambda_tilde)
}
sample_lambda.PG <- function(z, zit, fi, Time, timeidx, alpha, pri.lambda, W.sparse, W.dense){

  n <- length(fi)

  if(pri.lambda$type=="cps"|pri.lambda$type=="ind"){

    if(pri.lambda$type=="cps"){

      XWX <- sum(fi^2 * rowSums(W.dense)) # inner: summation over T; outer: summation over n
      inv.L0 <- 1/as.vector(pri.lambda$L0[1])
      Ln <- 1/(XWX+inv.L0)
      ln <- Ln*sum((zit*W.dense)*fi) # equivalent to Ln*sum(fi*rowSums(zitW))
      lambda <- rnorm(1, mean=ln, sd=sqrt(Ln))

    }else{ # independence prior

      W.tilde <- diag(colSums(fi^2*W.dense))
      Ln <- solve(W.tilde + diag(1/pri.lambda$L0, ncol(W.dense)))
      ln <- Ln %*% colSums(fi * (zit*W.dense))
      lambda <- MASS::mvrnorm(1, mu = ln, Sigma = Ln)

    }

  }

  if(pri.lambda$type=="rw1"|pri.lambda$type=="rw2"){

    # AWOL step

    zstar <- z-rep(alpha[1]*fi, Time)
    psi <- alpha[2]

    lambda_tilde <- AWOL_fac.PG(z = zstar, fi = fi, Time = Time,
                                psi = psi, pri.type = pri.lambda$type,
                                hyp.c = pri.lambda$c, W.dense = W.dense)

    lambda_tilde <- matrix(lambda_tilde, ncol = 1)

    # sample lambda and psi2

    A0lambda <- c(pri.lambda$phi, pri.lambda$zeta)

    Xf <- matrix(rep(fi, Time),ncol=1)
    if(pri.lambda$type=="rw1"){
      lth <- lambda_tilde[-1]
    } else{
      lth <- lambda_tilde
    }

    Zlambda <- constructZ(Xf, Time = Time, timeidx = timeidx, beta_tilde = lth)

    alpha_lambda <- sample_alpha.PG(z = z, d = 1, Z = Zlambda,
                                    A0 = A0lambda, W = W.sparse)

    lambdah <- transform_to_centered(beta_tilde = lambda_tilde,
                                     alpha = alpha_lambda,
                                     d = 1)
    if(pri.lambda$type=="rw1"){
      lambda<-lambdah[1+(1:Time)] # factor loading in observations
    }else{
      lambda<-lambdah
    }

    # sample phi_lambda and zeta_lambda

    phi_lambda <- sample_GIG(a = pri.lambda$a.phi,
                             l = pri.lambda$kappa.phi,
                             par = alpha_lambda[1])
    phi_lambda[phi_lambda>10^11]=10^11
    phi_lambda[phi_lambda<0.1^15]=0.1^15

    zeta_lambda <- sample_GIG(a = pri.lambda$a.zeta,
                              l = pri.lambda$kappa.zeta,
                              par = alpha_lambda[2])
    zeta_lambda[zeta_lambda>10^11]=10^11
    zeta_lambda[zeta_lambda<0.1^11]=0.1^11

    pri.lambda$phi <- phi_lambda
    pri.lambda$zeta <- zeta_lambda

    # sampling of kappa_phi and kappa_zeta

    if(pri.lambda$learn.kappa.phi){
      pri.lambda$kappa.phi <- sample_G(a = pri.lambda$a.phi,
                                       d = 1,
                                       par = pri.lambda$phi,
                                       prior_hp1 = pri.lambda$d.phi,
                                       prior_hp2 = pri.lambda$e.phi)
    }

    if(pri.lambda$learn.kappa.zeta){
      pri.lambda$kappa.zeta <- sample_G(a = pri.lambda$a.zeta,
                                        d = 1,
                                        par = pri.lambda$zeta,
                                        prior_hp1 = pri.lambda$d.zeta,
                                        prior_hp2 = pri.lambda$e.zeta)
    }

  }

  if(pri.lambda$type == "rw1" | pri.lambda$type == "rw2"){
    res <- list(lambda=lambda, alpha=alpha_lambda, pri.lambda=pri.lambda)
  } else{
    res <- list(lambda=lambda, alpha=alpha, pri.lambda=pri.lambda)
  }

  return(res)

}
AWOL_fac.PG <- function(z, fi, Time, psi, pri.type, hyp.c, W.dense){

  n <- length(fi)      # number of subjects
  Sf <- psi^2 * colSums(fi^2*W.dense)

  # construct prior matrix
  if(pri.type=="rw1"){
    d <- Time+1
    h <- list(rep(-1,d-1),
              c(1+1/hyp.c, rep(2, d-2),1))

    L0.inv <- as.matrix(Matrix::bandSparse(d, k=(-1):0, diag=h,symm=TRUE))
    Omega <- diag(c(0,Sf))+L0.inv
  }

  if(pri.type=="rw2"){
    d <- Time
    h <- list(rep(-1,d-1),
              c(1+1/hyp.c, rep(2, d-2),1))
    L0.inv <- as.matrix(Matrix::bandSparse(d, k=(-1):0, diag=h,symm=TRUE))
    Omega <- diag(Sf)+L0.inv
  }

  #create the c vector
  zh <- matrix(z, nrow = n, ncol = Time)
  cv <- matrix(psi*crossprod(fi, zh*W.dense), ncol = 1)

  if(pri.type=="rw1"){cv=c(0,cv)}

  # Oinv <- solve(Omega) added on 24.08.25
  Oinv <- MASS::ginv(Omega)
  m <- Oinv%*%cv
  lambda_tilde <- mvtnorm::rmvnorm(n=1, mean = m, sigma = Oinv)

  return(lambda_tilde)

}
