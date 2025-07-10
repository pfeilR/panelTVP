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
#' @details
#' Data for a panel with \eqn{n} subjects and Tmax \eqn{\equiv T} repeated measurements are
#'  simulated based on the linear predictor
#'  \deqn{\eta_{it} = \textbf{x}^\top_{it}\boldsymbol{\beta}_t + f_i \lambda_t, \quad i \in \{1,\dots,n\},
#'  \quad t \in \{1,\dots,T\},}
#'  where \eqn{\textbf{x}_{it}} is a column vector of the same dimension as the input vector
#'  \code{beta}. The first covariate is a 1 for estimating the global intercept, whereas
#'  the other covariates in \eqn{\textbf{x}} are independently generated standard Normal
#'  random variables. The subject-specific factor scores in \eqn{\textbf{f} = (f_1,\dots,f_n)^\top}
#'  are realizations of standard Normals as well. Moreover, factor scores are centered to have
#'  a mean of zero. The time-varying parameters
#'  \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and \eqn{\lambda_1,\dots,\lambda_T}
#'  are generated from first-order random walks. Using the non-centered parameterization
#'  (see Frühwirth-Schnatter and Wagner, 2010), the time-varying parameters are
#'  obtained by
#'  \deqn{\begin{aligned}
#'   \boldsymbol{\beta}_t &= \boldsymbol{\beta} + \boldsymbol{\Theta} \boldsymbol{\tilde{\beta}}_t, \\
#'   \lambda_t &= \lambda + \psi \tilde{\lambda}_t,
#'  \end{aligned}}
#'  where \eqn{\boldsymbol{\beta}} and \eqn{\lambda} are the fixed effects and correspond
#'  to the input arguments \code{beta} and \code{lambda}, respectively. The variation over
#'  time is controlled by \eqn{\boldsymbol{\Theta} = \text{diag}(\theta_1,\dots,\theta_d)} and
#'  \eqn{\psi}, which correspond to the input parameters \code{theta} and \code{psi},
#'  respectively. The parameters \eqn{\boldsymbol{\tilde{\beta}}_t, \tilde{\lambda}_t} are
#'  generated from standard Normal random walks.
#'
#'  Data for a \code{Gaussian} response are generated from
#'  \deqn{y_{it}|\eta_{it},\sigma^2 \sim \mathcal{N}(\eta_{it}, \sigma^2),}
#'  where \eqn{\sigma^2} is the homoscedastic error variance that is specified in the argument
#'  \code{sigma2}.
#'
#'  Data for a \code{Probit} response are generated from
#'  \deqn{y_{it}|\eta_{it} \sim \text{Bin}(1,\Phi(\eta_{it})),}
#'  where \eqn{\Phi(\cdot)} is the cumulative distribution function of the standard Normal distribution.
#'
#'  Data for a \code{Logit} response are generated from
#'  \deqn{y_{it}|\eta_{it}  \sim \text{Bin}(1,\frac{\exp(\eta_{it})}{1+\exp(\eta_{it})}).}
#'
#'  Data for a \code{Negative Binomial} response are generated from a Negative Binomial distribution
#'  that is parameterized (Pillow and Scott, 2012) as
#'  \deqn{p(y_{it}|r,\eta_{it}) = \frac{\Gamma(y_{it}+r)}{\Gamma(r)y_{it}!}(1-q_{it})^r q_{it}^{y_{it}},
#'   \quad q_{it} = \frac{\exp(\eta_{it})}{1+\exp(\eta_{it})},}
#'  where \eqn{r>0} is a common dispersion parameter that is specified in the argument \code{r}. Using
#'  this parameterization, the conditional moments are given by
#'  \deqn{\mathbb{E}(y_{it}|r,\eta_{it}) \equiv \mu_{it} = r \exp(\eta_{it}),\quad \mathbb{V}(y_{it}|r,\eta_{it}) = r \exp(\eta_{it}) (1+\exp(\eta_{it})).}
#'
#'  Data for a \code{Zero-Inflated Negative Binomial} (ZINB) response are generated according
#'  to the following mixture density (Neelon, 2019)
#'  \deqn{ y_{it}|r,\mu_{it},w_{it} \sim (1-\pi_{it}) \cdot \mathbb{I}_{(w_{it} = 0)} +
#'  \pi_{it} \cdot \mathcal{NB}(\mu_{it},r) \mathbb{I}_{(w_{it}=1)},}
#'  where \eqn{w_{it}} is a latent at-risk indicator, indicating whether or not a subject
#'  is at risk of experiencing the event. The at-risk indicators are generated as
#'  \deqn{w_{it}|\pi_{it} \sim \text{Bin}(1,\pi_{it}),}
#'  where the probability of at-risk membership is modelled via a Logit model as
#'  \deqn{\pi_{it} = \frac{\exp(\eta^\text{logit}_{it})}{1+\exp(\eta^\text{logit}_{it})}}
#'  Here, \eqn{\eta^\text{logit}_{it}} is a linear predictor with time-varying effects (defined as above).
#'  Observations that are at risk of experiencing the event are then simulated from a
#'  Negative Binomial model, where a separate linear predictor is used for modelling the count
#'  component of the ZINB model, i.e.,
#'  \deqn{y_{it}|r,\mu_{it},w_{it}=1 \sim \mathcal{NB}(\mu_{it}, r),}
#'  where \eqn{\mu_{it} = r \exp(\eta^\text{nb}_{it})}.
#'  Observations that are not at risk of experiencing the event, i.e., for which \eqn{w_{it}=0},
#'  are defined as structural zeros and, thus, set to zero.
#'  As the ZINB model contains two different linear predictors with time-varying parameters,
#'  simulation from this model requires the arguments \code{beta.nb}, \code{theta.nb}, \code{lambda.nb},
#'  \code{psi.nb}, \code{beta.logit}, \code{theta.logit}, \code{lambda.logit}, \code{psi.logit}.
#'  Those parameters serve the same purpose as the parameters \code{beta}, \code{theta}, \code{lambda},
#'  \code{psi} of the previous models. Furthermore, note that not the same sets of covariates have
#'  to be used in the zero-inflation and count predictors \eqn{\eta_{it}^\text{logit}} and
#'  \eqn{\eta_{it}^\text{nb}}, respectively. This simulation function generates different
#'  covariates for the two predictors, but the main function of this package [panelTVP()] can handle
#'  different sets.
#'
#' @description
#' This function simulates panel data with time-varying parameters,
#'  where the distribution of the
#'  response variable is either Gaussian, binary (using a Probit or Logit link)
#'  Negative Binomial or Zero-Inflated Negative Binomial.
#'
#' @references
#'  Frühwirth-Schnatter, S. and Wagner, H. (2010). Stochastic Model Specification
#'  Search for Gaussian and Partially Non-Gaussian State Space Models. Journal
#'  of Econometrics, 154, 85-100.
#'
#'  Knaus, P., Bitto-Nemling, A., Cadonna, A. and Frühwirth-Schnatter, S. (2021).
#'  Shrinkage in the Time-Varying Parameter Model Framework using the R Package
#'  \code{shrinkTVP}. In: Journal of Statistical Software, 100, 1-32.
#'
#'  Neelon, B. (2019). Bayesian Zero-Inflated Negative Binomial Regression Based
#'  on Pólya-Gamma Mixtures. In: Bayesian Analysis, 14, 829-855.
#'
#'  Pillow, J. and Scott, J. (2012). Fully Bayesian inference for neural models
#'  with negative-binomial spiking. In: Advances in neural information processing
#'  systems, 25.
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

