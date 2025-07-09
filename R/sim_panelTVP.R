#' Simulate data from a time-varying parameter panel data model
#'
#' @param n number of subjects (scalar)
#' @param Tmax number of time points / repeated measurements per subject (scalar)
#' @param model either "Gaussian", "Probit", "Logit", "NegBin" or "ZINB"
#' @param beta fixed regression effects of Gaussian, Probit, Logit or
#'  Negative Binomial model with the first value representing the
#'  global intercept (vector of dimension d)
#' @param theta standard deviations of random walk for regression effects
#'  of Gaussian, Probit, Logit or Negative Binomial model;
#'  larger values yield regression effects that vary stronger over time
#'  (vector of dimension d)
#' @param lambda fixed factor loading of Gaussian, Probit, Logit or
#'  Negative Binomial model (scalar)
#' @param psi standard deviation of random walk for factor loading
#'  of Gaussian, Probit, Logit or Negative Binomial model;
#'  a larger value yields a factor loading that varies stronger over time (scalar)
#' @param beta.nb fixed regression effects in count component of Zero-Inflated
#'  Negative Binomial model with the first value representing the
#'  global intercept (vector of dimension d_nb)
#' @param theta.nb standard deviations of random walk for regression effects
#'  in count component of Zero-Inflated Negative Binomial model;
#'  larger values yield regression effects that vary stronger over time
#'  (vector of dimension d_nb)
#' @param lambda.nb fixed factor loading in count component
#'  of Zero-Inflated Negative Binomial model (scalar)
#' @param psi.nb standard deviation of random walk for factor loading
#'  in count component Zero-Inflated Negative Binomial model;
#'  a larger value yields a factor loading that varies stronger over time (scalar)
#' @param beta.logit fixed regression effects in zero-inflation component of Zero-Inflated
#'  Negative Binomial model with the first value representing the
#'  global intercept (vector of dimension d_logit)
#' @param theta.logit standard deviations of random walk for regression effects
#'  in zero-inflation component of Zero-Inflated Negative Binomial model;
#'  larger values yield regression effects that vary stronger over time
#'  (vector of dimension d_logit)
#' @param lambda.logit fixed factor loading in zero-inflation component
#'  of Zero-Inflated Negative Binomial model (scalar)
#' @param psi.logit standard deviation of random walk for factor loading
#'  in zero-inflation component Zero-Inflated Negative Binomial model;
#'  a larger value yields a factor loading that varies stronger over time (scalar)
#' @param r dispersion parameter of the Negative Binomial and
#'  Zero-Inflated Negative Binomial model (scalar, ignored
#'  if response is not simulated from a Negative Binomial or
#'  Zero-Inflated Negative Binomial model)
#' @param sigma2 homoscedastic error variance of the Gaussian model (scalar,
#'  ignored if response is not simulated from a Gaussian regression model)
#'
#' @returns a list of simulated data and parameters that contains the following elements:
#' \describe{
#'  \item{observed}{a data frame that contains the response variable
#'   \code{y}, the covariates (starting with the letter \code{W}),
#'   the time index \code{t} and the subject index \code{id}; for the ZINB model,
#'   the covariates of the zero-inflation component start with the letter \code{W}
#'   and end with \code{logit}, whereas the covariates of the count component
#'   start with the letter \code{W} and end with \code{nb})}
#'  \item{beta}{a T x d matrix of regression effects, i.e., each row contains
#'   the regression effects of the corresponding time point (only in Gaussian,
#'   Probit, Logit and Negative Binomial model)}
#'  \item{lambda}{a T x 1 matrix of factor loadings, i.e., each row contains
#'   the factor loading of the corresponding time point (only in Gaussian,
#'   Probit, Logit and Negative Binomial model)}
#'   \item{beta.logit}{a T x d_logit matrix of regression effects
#'    in zero-inflation component of the model, i.e., each row contains
#'    the regression effects of the corresponding time point (only in ZINB model)}
#'   \item{beta.nb}{a T x d_nb matrix of regression effects
#'    in count component of the model, i.e., each row contains
#'    the regression effects of the corresponding time point (only in ZINB model)}
#'  \item{lambda.logit}{a T x 1 matrix of factor loadings in
#'   zero-inflation component of the model, i.e., each row contains
#'   the factor loading of the corresponding time point (only in ZINB model)}
#'   \item{lambda.nb}{a T x 1 matrix of factor loadings in
#'    count component of the model, i.e., each row contains
#'    the factor loading of the corresponding time point (only in ZINB model)}
#'   \item{sigma2}{a scalar containing the value of \eqn{\sigma^2} (only in
#'   Gaussian model)}
#'   \item{r}{a scalar containing the value of \eqn{r} (only in Negative
#'   Binomial and ZINB model)}
#' }
#'
#' @description
#' This function simulates panel data with time-varying parameters,
#'  where the distribution of the
#'  response variable is either Gaussian, binary (using a Probit or Logit link)
#'  Negative Binomial or Zero-Inflated Negative Binomial.
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
#'
#' # Simulating data from a Zero-Inflated Negative Binomial panel model
#' x <- sim_panelTVP(n = 100, Tmax = 4,
#'                   beta.nb = c(0.5,0,-0.2), theta.nb = c(0.1,0,0),
#'                   lambda.nb = 0.3, psi.nb = 0,
#'                   beta.logit = c(0.9,0.1,0), theta.logit = c(0.05,0,0),
#'                   lambda.logit = 0.8, psi.logit = 0.1,
#'                   model = "ZINB", r = 2)
#' head(x$observed, 10)
#' x$beta.logit
#' x$beta.nb
#' x$lambda.logit
#' x$lambda.nb
#' x$r
#' @author Roman Pfeiler, Helga Wagner
#' @export
sim_panelTVP <- function(n,
                         Tmax,
                         model,
                         beta = NULL,
                         theta = NULL,
                         lambda = NULL,
                         psi = NULL,
                         beta.nb = NULL,
                         theta.nb = NULL,
                         lambda.nb = NULL,
                         psi.nb = NULL,
                         beta.logit = NULL,
                         theta.logit = NULL,
                         lambda.logit = NULL,
                         psi.logit = NULL,
                         r = NULL,
                         sigma2 = NULL){

  # Input Checks ---------------------------------------------------------------

  check_sim(n = n,
            Tmax = Tmax,
            model = model,
            beta = beta,
            theta = theta,
            lambda = lambda,
            psi = psi,
            r = r,
            sigma2 = sigma2,
            beta.nb = beta.nb,
            theta.nb = theta.nb,
            lambda.nb = lambda.nb,
            psi.nb = psi.nb,
            beta.logit = beta.logit,
            theta.logit = theta.logit,
            lambda.logit = lambda.logit,
            psi.logit = psi.logit)

  # Simulate Data --------------------------------------------------------------

  if(!(model %in% "ZINB")){
    result <- sim_Gaussian_Probit_Logit_NegBin(n = n,
                                               Tmax = Tmax,
                                               beta = beta,
                                               theta = theta,
                                               lambda = lambda,
                                               psi = psi,
                                               model = model,
                                               r = r,
                                               sigma2 = sigma2)

  } else{
    result <- sim_ZINB(n = n,
                       Tmax = Tmax,
                       beta.nb = beta.nb,
                       theta.nb = theta.nb,
                       lambda.nb = lambda.nb,
                       psi.nb = psi.nb,
                       beta.logit = beta.logit,
                       theta.logit = theta.logit,
                       lambda.logit = lambda.logit,
                       psi.logit = psi.logit,
                       r = r)
  }

  return(result)

}

