#' Simulate data from a time-varying parameter panel data model
#'
#' @param n number of subjects (scalar)
#' @param Tmax number of time points / repeated measurements per subject (scalar)
#' @param beta fixed regression effects, with the first value representing the
#'  global intercept (vector of dimension d)
#' @param theta standard deviation of random walk for regression effects, i.e.,
#'  larger values yield regression effects that vary stronger over time
#'  (vector of dimension d)
#' @param lambda fixed factor loading (scalar)
#' @param psi standard deviation of random walk for factor loading, i.e.,
#'  a larger value yields a factor loading that varies stronger over time (scalar)
#' @param model either "Gaussian", "Probit", "Logit" or "NegBin"
#' @param r dispersion parameter of the Negative Binomial model (scalar, ignored
#'  if response is not simulated from a Negative Binomial regression model)
#' @param sigma2 homoscedastic error variance of the Gaussian model (scalar,
#'  ignored if response is not simulated from a Gaussian regression model)
#'
#' @returns a list of simulated data and parameters that contains the following elements:
#' \describe{
#'  \item{observed}{a data frame that contains the response variable
#'   \code{y} the covariates (starting with the letter \code{W}),
#'   the time index \code{t} and the subject index \code{id}}
#'  \item{beta}{a T x d matrix of regression effects, i.e., each row contains
#'   the regression effects of the corresponding time point}
#'  \item{lambda}{a T x 1 matrix of factor loadings, i.e., each row contains
#'   the factor loading of the corresponding time point}
#'   \item{sigma2}{a scalar containing the value of \eqn{\sigma^2} (only in
#'   Gaussian model)}
#'   \item{r}{a scalar containing the value of \eqn{r} (only in Negative
#'   Binomial model)}
#' }
#'
#' @description
#' This function simulates panel data with time-varying parameters,
#'  where the distribution of the
#'  response variable is either Gaussian, binary (using a Probit or Logit link)
#'  or Negative Binomial.
#'
#' @examples
#' # Simulating data from a Gaussian panel model
#' x <- sim_panelTVP(n = 100, Tmax = 6,
#'                   beta = c(4,1,0), theta = c(0,1,0),
#'                   lambda = 1, psi = 0.2,
#'                   model = "Gaussian", sigma2 = 1)
#' head(x$observed, 10)
#' x$beta
#' x$lambda
#' x$sigma2
#'
#' # Simulating data from a Probit panel model
#' x <- sim_panelTVP(n = 100, Tmax = 6,
#'                   beta = c(0.5,0.2,0), theta = c(0.1,0.25,0),
#'                   lambda = 0.1, psi = 0,
#'                   model = "Probit")
#' head(x$observed, 10)
#' x$beta
#' x$lambda
#'
#' # Simulating data from a Logit panel model
#' x <- sim_panelTVP(n = 100, Tmax = 6,
#'                   beta = c(0.5,0.2,0), theta = c(0.1,0.25,0),
#'                   lambda = 0.5, psi = 0.04,
#'                   model = "Logit")
#' head(x$observed, 10)
#' x$beta
#' x$lambda
#'
#' # Simulating data from a Negative Binomial panel model
#' x <- sim_panelTVP(n = 100, Tmax = 6,
#'                   beta = c(0.5,0.2,0), theta = c(0.1,0.25,0),
#'                   lambda = 0.9, psi = 0.1,
#'                   model = "NegBin", r = 2)
#' head(x$observed, 10)
#' x$beta
#' x$lambda
#' x$r
#' @export
sim_panelTVP <- function(n,
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
  Sigma <- diag(rep(sigma2, t))
  if(model == "Gaussian") error <- matrix(mvtnorm::rmvnorm(n, mean = rep(0, t), sigma = Sigma))
  y <- matrix(NA,n*t)
  fi <- matrix(rnorm(n), ncol = 1)
  fi <- fi - mean(fi)
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

  if(model == "Gaussian"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                sigma2 = sigma2,
                model = "Gaussian")

  }

  if(model == "Probit"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                model = "Probit")

  }

  if(model == "Logit"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                model = "Logit")

  }

  if(model == "NegBin"){

    ret <- list(observed = observed,
                beta = betat,
                lambda = lambdat,
                r = r,
                model = "Negative Binomial")

  }

  return(ret)

}
