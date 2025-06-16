#' Simulate data from a time-varying parameter panel data model
#'
#' @param n number of subjects (scalar)
#' @param t number of time points (scalar)
#' @param beta starting values of random walk prior for regression effects
#'  (vector of dimension d)
#' @param theta standard deviation of random walk prior for regression effects
#'  (vector of dimension d)
#' @param lambda starting value of random walk prior for factor loading (scalar)
#' @param psi standard deviation of random walk prior for factor loading (scalar)
#' @param model either "Gaussian", "Probit", "Logit" or "NegBin"
#' @param r dispersion parameter of the Negative Binomial model (scalar, ignored
#'  if response is not simulated from a Negative Binomial regression model)
#' @param sigma2 homoscedastic error variance of the Gaussian model (scalar,
#'  ignored if response is not simulated from a Gaussian regression model)
#'
#' @returns a list with simulated data and parameters
#' @export
#'
#' @examples 1
#'
#' @description
#' This function simulates balanced panel data, where the distribution of the
#'  response variable is either Gaussian, binary (using a Probit or Logit link)
#'  or Negative Binomial.
#'
sim_panelTVP <- function(n,
                         t,
                         beta = NULL,
                         theta = NULL,
                         lambda = NULL,
                         psi = NULL,
                         model,
                         r = NULL,
                         sigma2 = NULL){

  # input checks -> To-Do!

  # if(!is.numeric(n) | length(n)>1){
  #   stop("n needs to be a numeric of length 1")
  # }
  # if(!is.numeric(t) | length(t)>1){
  #   stop("t needs to be a numeric of length 1")
  # }
  # if(is.null(beta) & is.null(betat)){
  #   stop("you need to specify either (beta and theta) or betat")
  # }
  # if(is.null(theta) & is.null(betat)){
  #   stop("you need to specify either (beta and theta) or betat")
  # }
  # if(!is.vector(beta) | !is.vector(theta)){
  #   stop("beta and theta have to be of vectors")
  # }
  # if(!is.null(betat)){
  #   if(!is.matrix(betat)){
  #     stop("betat has to be a matrix of dimension t x d")
  #   }
  #   if(nrow(betat)!=t){
  #     stop("betat has to be a matrix of dimension t x d")
  #   }
  # }
  # if(is.null(lambda) & is.null(lambdat)){
  #   stop("you need to specify either (lambda and psi) or lambdat")
  # }
  # if(is.null(psi) & is.null(lambdat)){
  #   stop("you need to specify either (lambda and psi) or lambdat")
  # }
  # if(!is.vector(lambda) | !is.vector(psi)){
  #   stop("lambda and psi have to be numerics each of length 1")
  # }
  # if(!is.null(lambdat)){
  #   if(!is.matrix(lambdat)){
  #     stop("lambdat has to be a matrix of dimension t x 1")
  #   }
  #   if(nrow(lambdat)!=t){
  #     stop("lambdat has to be a matrix of dimension t x 1")
  #   }
  #   if(ncol(lambdat)!=1){
  #     stop("lambdat has to be a matrix of dimension t x 1")
  #   }
  # }
  # if(!(model %in% c("Gaussian", "Probit", "Logit", "NegBin"))){
  #   stop("model needs to be either Gaussian, Probit, Logit or NegBin - no default!")
  # }
  # if(model == "NegBin" & is.null(r)){
  #   stop("r needs to be specified when simulating from a Negative Binomial model")
  # }
  # if(!is.null(sigma2) & (!is.numeric(sigma2) | length(sigma2)>1)){
  #   stop("sigma2 needs to be a numeric of length 1")
  # }
  # if(model != "Gaussian" & !is.null(sigma2)){
  #   warning("sigma2 is ignored as data from a non-Gaussian model are simulated")
  # }
  # if(model == "Gaussian" & is.null(sigma2)){
  #   stop("sigma2 needs to be specified when simulating from a Gaussian model")
  # }
  # if(!is.null(betat)){
  #   if(ncol(betat) < 2) stop("we need at least 2 covariates (intercept + feature)")
  # }
  # if(!is.null(beta)){
  #   if(length(beta) != length(theta)){
  #     stop("beta and theta need to be of the same length")
  #   }
  #   if(length(beta) < 2) stop("we need at least 2 covariates (intercept + feature)")
  # }
  # if(!is.null(lambda)){
  #   if(length(lambda) != 1) stop("lambda has to be a scalar")
  #   if(length(psi) != 1) stop("psi has to be a scalar")
  # }

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
  Sigma <- diag(rep(sigma2, t))
  if(model == "Gaussian") error <- matrix(mvtnorm::rmvnorm(n, mean = rep(0, t), sigma = Sigma))
  y <- matrix(NA,n*t)
  fi <- matrix(rnorm(n), ncol = 1)
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

    if(model == "Probit"){
      prob.binary <- pnorm(eta)
      y <- rbinom(length(prob.binary), size = 1, prob = prob.binary)
    }

    if(model == "NegBin"){
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


  if(model == "NegBin"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                r = r)

  } else{

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat)

  }

  return(ret)

}
