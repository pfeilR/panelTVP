#' Simulate data from a Zero-Inflated Negative Binomial (ZINB) panel with
#'  time-varying parameters.
#'
#' @param n number of subjects (scalar)
#' @param Tmax number of time points / repeated measurements per subject (scalar)
#' @param beta.nb fixed regression effects, with the first value representing the
#'  global intercept (vector of dimension d_nb) (count model)
#' @param theta.nb standard deviation of random walk for regression effects, i.e.,
#'  larger values yield regression effects that vary stronger over time
#'  (vector of dimension d_nb) (count model)
#' @param lambda.nb fixed factor loading (scalar) (count model)
#' @param psi.nb standard deviation of random walk for factor loading, i.e.,
#'  a larger value yields a factor loading that varies stronger over time (scalar)
#'  (count model)
#' @param beta.logit fixed regression effects, with the first value representing the
#'  global intercept (vector of dimension d_logit) (zero-inflation model)
#' @param theta.logit standard deviation of random walk for regression effects, i.e.,
#'  larger values yield regression effects that vary stronger over time
#'  (vector of dimension d_logit) (zero-inflation model)
#' @param lambda.logit fixed factor loading (scalar) (zero-inflation model)
#' @param psi.logit standard deviation of random walk for factor loading, i.e.,
#'  a larger value yields a factor loading that varies stronger over time (scalar)
#'  (zero-inflation model)
#' @param r dispersion parameter of the Negative Binomial component of the ZINB model
#'
#' @returns a list of simulated data and parameters that contains the following elements:
#' \describe{
#'  \item{observed}{a data frame that contains the response variable
#'   \code{y} the covariates of the zero-inflation component
#'    (starting with the letter \code{W} and ending with \code{logit}), the
#'    covariates of the count component (starting with the letter \code{W} and
#'    ending with \code{nb}), the time index \code{t} and the subject index \code{id}}
#'  \item{beta.logit}{a T x d_logit matrix of regression effects
#'    for the zero-inflation component of the model, i.e., each row contains
#'    the regression effects of the corresponding time point}
#'   \item{beta.nb}{a T x d_nb matrix of regression effects
#'    for the count component of the model, i.e., each row contains
#'    the regression effects of the corresponding time point}
#'  \item{lambda.logit}{a T x 1 matrix of factor loadings for the
#'   zero-inflation component of the model, i.e., each row contains
#'   the factor loading of the corresponding time point}
#'   \item{lambda.nb}{a T x 1 matrix of factor loadings for the
#'    count component of the model, i.e., each row contains
#'    the factor loading of the corresponding time point}
#'   \item{r}{a scalar containing the value of \eqn{r} for the count component of the model}
#' }
#' @export
#'
#' @examples
#' x <- sim_panelTVP_ZINB(n = 100,
#'                        Tmax = 4,
#'                        beta.nb = c(0.5,0,-0.2),
#'                        theta.nb = c(0.1,0,0),
#'                        lambda.nb = 0.3,
#'                        psi.nb = 0,
#'                        beta.logit = c(0.9,0.1,0),
#'                        theta.logit = c(0.05,0,0),
#'                        lambda.logit = 0.8,
#'                        psi.logit = 0.1,
#'                        r = 2)
#' head(x$observed, 10)
#' x$beta.logit
#' x$beta.nb
#' x$lambda.logit
#' x$lambda.nb
#' x$r
sim_panelTVP_ZINB <- function(n,
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
  fi.logit <- matrix(rnorm(n), ncol = 1)
  fi.logit <- fi.logit - mean(fi.logit)
  Fmat.logit <- matrix(NA,n*t)
  fi.nb <- matrix(rnorm(n), ncol = 1)
  fi.nb <- fi.nb - mean(fi.nb)
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
              beta.logit = betat.logit,
              beta.nb = betat.nb,
              lambda.logit = lambdat.logit,
              lambda.nb = lambdat.nb,
              r = r)
  return(ret)

}
