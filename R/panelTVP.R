#' Fit a Bayesian panel model with time-varying parameters
#'
#' @param formula the usual formula argument in regression methods, e.g., as in [lm()]
#'  (no default)
#' @param data a data.frame containing the variables specified in the formula argument
#'  (no default)
#' @param model a character indicating which model you want to estimate
#'  This parameter is either 'Gaussian', 'Probit', 'Logit' or 'NegBin' depending on
#'  whether you want to fit a model for Gaussian, Probit, Logit or Negative Binomial
#'  response data. In case you want to estimate a Zero-Inflated Negative Binomial
#'  model, please use the function [`panelTVP_ZINB()`], which was designed for
#'  this purpose
#' @param prior.reg a list of arguments for estimating the parameters of the regression
#'  part of the model. The arguments are:
#'   \itemize{
#'    \item \code{d.tau}: shape parameter of Gamma prior for \eqn{\kappa^\tau}
#'    \item \code{e.tau}: rate parameter of Gamma prior for \eqn{\kappa^\tau}
#'    \item \code{d.xi}: shape parameter of Gamma prior for \eqn{\kappa^\xi}
#'    \item \code{e.xi}: rate parameter of Gamma prior for \eqn{\kappa^\xi}
#'    \item \code{b.tau}: part of the rate parameter of the Gamma prior for \eqn{a^\tau}
#'    \item \code{nu.tau}: shape parameter of the Gamma prior for \eqn{a^\tau}
#'    \item \code{b.xi}: part of the rate parameter of the Gamma prior for \eqn{a^\xi}
#'    \item \code{nu.xi}: shape parameter of the Gamma prior for \eqn{a^\xi}
#'    \item \code{a.tau}: shape parameter of the Gamma prior for \eqn{\tau^2_j}
#'    \item \code{kappa.tau}: part of the rate parameter of the Gamma prior for \eqn{\tau^2_j}
#'    \item \code{a.xi}: shape parameter of the Gamma prior for \eqn{\xi^2_j}
#'    \item \code{kappa.xi}: part of the rate parameter of the Gamma prior for \eqn{\xi^2_j}
#'    \item \code{iota.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\tau}
#'    \item \code{iota.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\xi}
#'    \item \code{learn.a.tau}: if TRUE \eqn{a^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.tau} used as starting value, if FALSE
#'     \eqn{a^\tau =} \code{a.tau}
#'    \item \code{learn.a.xi}: if TRUE \eqn{a^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.xi} used as starting value, if FALSE
#'     \eqn{a^\xi =} \code{a.xi}
#'    \item \code{tau.target.rate}: desired acceptance rate when updating
#'     \eqn{a^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.tau = FALSE})
#'    \item \code{xi.target.rate}: desired acceptance rate when updating
#'     \eqn{a^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.xi = FALSE}
#'    \item \code{learn.kappa.tau}: if TRUE \eqn{\kappa^\tau} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.tau} used as starting value,
#'     if FALSE \eqn{\kappa^\tau = } \code{kappa.tau}
#'    \item \code{learn.kappa.xi}: if TRUE \eqn{\kappa^\xi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.xi} used as starting value,
#'     if FALSE \eqn{\kappa^\xi = } \code{kappa.xi}
#'    \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "rw1" (RW1 shrinkage prior), "rw2" (RW2 shrinkage prior)
#'      or "ind" (independence prior)
#'    \item \code{c}: prior parameter that scales the variance when using the
#'     random walk shrinkage prior (ignored when using independence prior)
#'    \item \code{B0}: prior variance on the regression parameters when using the
#'     independence prior (ignored when using shrinkage prior)
#'   }
#' @param prior.var a list of arguments for estimating the homoscedastic error variance
#'  in a Gaussian/Normal model. For other models, this argument is ignored.
#'  The arguments are:
#'  \itemize{
#'   \item \code{learn.C0.hyp}: this argument is a list containing the prior
#'    parameters for the prior variance \eqn{C_0} of \eqn{\sigma^2}
#'   \item \code{c0}: variance parameter of the Inverse-Gamma prior on \eqn{\sigma^2}
#'    }
#' @param prior.load a list of arguments for estimating the parameters of the factor
#'  part of the model. The arguments are:
#'  \itemize{
#'   \item \code{d.phi}: shape parameter of Gamma prior for \eqn{\kappa^\phi}
#'   \item \code{e.phi}: rate parameter of Gamma prior for \eqn{\kappa^\phi}
#'   \item \code{d.zeta}: shape parameter of Gamma prior for \eqn{\kappa^\zeta}
#'   \item \code{e.zeta}: rate parameter of Gamma prior for \eqn{\kappa^\zeta}
#'   \item \code{a.phi}: shape parameter of the Gamma prior for \eqn{\phi^2_j}
#'   \item \code{kappa.phi}: part of the rate parameter of the Gamma prior for \eqn{\phi^2_j}
#'   \item \code{a.zeta}: shape parameter of the Gamma prior for \eqn{\zeta^2_j}
#'   \item \code{kappa.zeta}: part of the rate parameter of the Gamma prior for \eqn{\zeta^2_j}
#'   \item \code{learn.kappa.phi}: if TRUE \eqn{\kappa^\phi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.phi} used as starting value,
#'     if FALSE \eqn{\kappa^\phi = } \code{kappa.phi}
#'   \item \code{learn.kappa_zeta}: if TRUE \eqn{\kappa^\zeta} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.zeta} used as starting value,
#'     if FALSE \eqn{\kappa^\zeta = } \code{kappa.zeta}
#'   \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "cps" (compound symmetric), "rw1" (RW1 shrinkage prior),
#'     "rw2" (RW2 shrinkage prior) or "ind" (independence prior)
#'   \item \code{c}: prior parameter that scales the variance when using the
#'     random walk shrinkage prior (ignored when using cps or independence prior)
#'   \item \code{L0}: prior variance on the factor loading (ignored when using shrinkage prior)
#'  }
#'  Note that for the factor model, the hyperparameters \code{a_phi} and
#'  \code{a_zeta} have to be fixed and are not sampled using Metropolis-Hastings
#' @param mcmc.opt a list containing information for the overall sampler.
#'  The arguments are:
#'  \itemize{
#'   \item \code{chain.length}: the length of the Markov Chain (i.e., total number
#'   of draws)
#'   \item \code{burnin}: the burn-in period
#'   \item \code{thin}: the thinning factor
#'   \item \code{asis}: if set to TRUE, an ancillarity sufficiency interweaving
#'    step is added for increasing the sampling efficiency of the regression
#'    effects
#'  }
#'  Note that the final Markov Chain (after applying burn-in and thinning)
#'  is of length \deqn{\frac{\text{chain.length}-\text{burnin}}{\text{thin}}}
#' @param settings.NegBin a list containing information for sampling the dispersion
#'  parameter r in the Negative Binomial model using univariate Slice sampling.
#'  For other response distributions, this is ignored. The arguments are:
#'  \itemize{
#'    \item \code{alpha.r}: shape parameter of Gamma proposal
#'    \item \code{beta.r}: rate parameter of Gamma proposal
#'    \item \code{expansion.steps}: number of steps in stepping-out phase
#'    \item \code{width}: width of the slice interval
#'    \item \code{p.overrelax}: probability of performing an overrelaxation step;
#'     overrelaxation may help boosting sampling efficiency; when overrelaxation
#'     should not be used, set \code{p.overrelax = 0}
#'    \item \code{accuracy.overrelax}: accuracy in overrelaxation phase
#'  }
#'  For more information on overrelaxation and Slice sampling in general,
#'  see the original paper on Slice sampling by Neal (2003)
#' @param HPD.coverage coverage probability of highest posterior density intervals
#'  (default yields 95 percent coverage)
#'
#' @description
#' This is the main function for fitting a flexible Bayesian panel data model
#'  in the time-varying parameter framework. By using shrinkage priors, it is
#'  possible to identify whether an effect is time-varying, time-invariant or zero.
#'  This function works for Gaussian, binary and Negative Binomial response data.
#'  When dealing with count data subject to zero-inflation, the function
#'  [panelTVP_ZINB()] should be used instead.
#'
#' @details
#'  This function fits a Bayesian time-varying parameter panel data model to a
#'  longitudinal response \eqn{y_{it}} for \eqn{i \in \{1,\dots,n\}} subjects that are observed
#'  at \eqn{t \in \{1,\dots,\text{T}\}}
#'  time point. The model is expressed in its non-centered form
#'  (see Frühwirth-Schnatter and Wagner (2010) for details on this parameterization)
#'  as this makes it easier to place shrinkage priors on the model parameters.
#'  By using the non-centered parameterization, the response is linked to the
#'  following linear predictor
#'  \deqn{\eta_{it} = \mathbf{x}_{it}^\top \boldsymbol{\beta} +
#'  \mathbf{x}_{it}^\top \boldsymbol{\Theta} \boldsymbol{\tilde{\beta}}_t +
#'   f_i\lambda + f_i \psi \tilde{\lambda}_t,}
#'  where \eqn{\mathbf{x}_{it}} is a d-dimensional column vector of covariate
#'  values for subject \eqn{i} at time \eqn{t}, \eqn{\boldsymbol{\beta}} is a d-dimensional
#'  column vector of fixed effects (including the intercept as first parameter),
#'  \eqn{\boldsymbol{\Theta} = \text{diag}(\theta_1,\dots,\theta_d)} is a
#'  diagonal matrix, where larger main diagonal elements indicate stronger variation
#'  of the regression effects over time, \eqn{\boldsymbol{\tilde{\beta}}_t}
#'  is a d-dimensional state vector that follows a Gaussian random walk,
#'  \eqn{f_i} is a subject-specific factor score, \eqn{\lambda} is the fixed loading,
#'   \eqn{\psi} is the standard deviation of the factor loading, where a larger value
#'   indicates stronger variation of the factor loading over time, and
#'   \eqn{\tilde{\lambda}_t} is a scalar that follows a standard Normal random walk.
#'
#'  The model only requires priors on the parameters
#'  \eqn{\boldsymbol{\beta}, \lambda, \psi} and the diagonal elements of
#'  \eqn{\boldsymbol{\Theta}}. By using the non-centered parameterization of the
#'  model, standard priors from regression analysis can be used for all parameters.
#'  Following Bitto and Frühwirth-Schnatter (2019), we place hierarchical
#'  Normal-Gamma shrinkage priors on those parameters. This makes it possible to
#'  identify whether an effect is time-varying, time-invariant or zero and, thus,
#'  prevents the model from overfitting (i.e., it is reasonable to assume that not
#'  every covariate has a time-varying effect and without proper regularization
#'  estimates are likely less stable).
#'  The priors on the parameters are specified as
#'  \deqn{
#'  \begin{aligned}
#'    \theta_j|\xi_j^2 &\sim \mathcal{N}(0,\xi^2_j),
#'    &\quad \xi_j^2|a^\xi, \kappa^\xi &\sim \mathcal{G}\left(a^\xi, \frac{a^\xi \kappa^\xi}{2}\right),
#'    &\quad j = \{1, \dots, d\}, \\
#'    \beta_j|\tau_j^2 &\sim \mathcal{N}(0,\tau_j^2),
#'    &\quad \tau_j^2|a^\tau, \kappa^\tau &\sim \mathcal{G}\left(a^\tau, \frac{a^\tau \kappa^\tau}{2}\right),
#'    &\quad j = \{1, \dots, d\}, \\
#'    \psi|\zeta^2 &\sim \mathcal{N}(0,\zeta^2),
#'    &\quad \zeta^2|a^\zeta, \kappa^\zeta &\sim \mathcal{G}\left(a^\zeta, \frac{a^\zeta \kappa^\zeta}{2}\right), \\
#'    \lambda|\phi^2 &\sim \mathcal{N}(0,\phi^2),
#'    &\quad \phi^2|a^\phi, \kappa^\phi &\sim \mathcal{G}\left(a^\phi, \frac{a^\phi \kappa^\phi}{2}\right).
#'  \end{aligned}
#'  }
#'  The hyperparameters \eqn{a^\zeta,a^\phi} in the factor part of the model are
#'  held at fixed values (as in our simulations we have experienced instability issues
#'  when sampling those parameters of the factor model as well),
#'  whereas all the other hyperparameters may either be held
#'  fixed or may be equipped with additional hyperpriors. In the latter case the
#'  following hyperpriors are considered:
#'  \deqn{
#'  \begin{aligned}
#'   \kappa^\xi|d^\xi,e^\xi& \sim \mathcal{G}(d^\xi,e^\xi), \quad  a^\xi|\nu^\xi,b^\xi \sim  \mathcal{G}(\nu^\xi, \nu^\xi b^\xi), \\
#'   \kappa^\tau|d^\tau,e^\tau& \sim  \mathcal{G}(d^\tau,e^\tau), \quad a^\tau|\nu^\tau,b^\tau \sim  \mathcal{G}(\nu^\tau, \nu^\tau b^\tau), \\
#'   \kappa^\zeta|d^\zeta,e^\zeta& \sim \mathcal{G}(d^\zeta,e^\zeta), \\
#'   \kappa^\phi|d^\phi,e^\phi& \sim \mathcal{G}(d^\phi,e^\phi).
#'  \end{aligned}
#' }
#'  Note that for \eqn{a \le 1} the priors are valid shrinkage priors. The
#'  Bayesian Lasso is a special case of the Normal-Gamma prior with \eqn{a = 1}.
#'  In our simulations we have achieved good results by using the Bayesian Lasso
#'  for the hyperparameters in the factor model (\eqn{a^\phi, a^\zeta}) and by sampling
#'  the hyperparameters in the regression model (\eqn{a^\tau, a^\xi}) using
#'  Metropolis-Hastings updates. However, we recommend to carry out a
#'  prior sensitivity analysis for the fixed hyperparameters.
#'
#'  The function \code{panelTVP} can handle the following popular regression models
#'   \itemize{
#'     \item Gaussian model (Normal outcomes)
#'     \item Logit model (binary outcomes)
#'     \item Probit model (binary outcomes)
#'     \item Negative Binomial model (overdispersed count outcomes)
#'   }
#'
#'  With a Gaussian (Normal) response, the model is given as follows
#'  \deqn{y_{it} = \eta_{it} + \varepsilon_{it}, \quad \varepsilon_{it} \sim \mathcal{N}(0,\sigma^2),}
#'  where we place a hierarchical prior on the error variance
#'  \deqn{\sigma^2|c_0,C_0 \sim \mathcal{G}^{-1}(c_0,C_0), \quad C_0|g_0,G_0 \sim \mathcal{G}(g_0,G_0).}
#'
#'  In the binary Logit model, the probability that \eqn{y_{it} = 1} is modelled
#'  via the Logit-link as
#'  \deqn{\mathbb{P}(y_{it} = 1) = \frac{\exp(\eta_{it})}{1+\exp(\eta_{it})},   }
#'  whereas in the binary Probit model, this probability is modelled
#'  via the Probit-link as
#'  \deqn{\mathbb{P}(y_{it} = 1) = \Phi(\eta_{it}),}
#'  where \eqn{\Phi(\cdot)} denotes the standard Normal cumulative distribution function.
#'
#'  For modelling overdispersed count data, we assume that the response is a realization
#'  of a Negative Binomial distribution with pmf given by
#'  \deqn{p(y_{it}|r,\eta_{it}) = \frac{\Gamma(y_{it}+r)}{\Gamma(r)y_{it}!}(1-q_{it})^r q_{it}^{y_{it}},
#'   \quad q_{it} = \frac{\exp(\eta_{it})}{1+\exp{\eta_{it}}},}
#'  where we place a prior on the constant dispersion parameter
#'   \deqn{r|\alpha^r,\beta^r \sim \mathcal{G}(\alpha^r,\beta^r).}
#'  Note that we use this parameterization of the Negative Binomial distribution
#'  to facilitate posterior inference with Pólya-Gamma random variables
#'  (see Pillow and Scott, 2012).
#'
#' @returns The function returns an object of class \code{"panelTVP.Gaussian"},
#'  \code{"panelTVP.Probit"}, \code{"panelTVP.Logit"} or \code{"panelTVP.NegBin"}
#'  depending on whether a Gaussian, Probit, Logit or Negative Binomial model was
#'  fitted. The returned object contains a list of the following elements:
#'  \describe{
#'    \item{data}{the data used for fitting the model and additional context information
#'    derived from the data}
#'    \item{Y}{the \eqn{Tn \times M} response data matrix of every iteration of the chain,
#'    i.e., this matrix is only included in the output when missing response data
#'    were present and imputed via data augmentation. Each row of \code{Y} contains
#'    the sampled values of a specific observation. For observed data, the corresponding
#'    rows are essentially replicates of the same value. For missing data, the
#'    corresponding rows contain the imputed values of every iteration.}
#'    \item{mcmc}{Markov Chains for every parameter except for the factor scores
#'     (to save memory)}
#'    \item{posterior}{preliminary summary of posterior results}
#'    \item{fmean}{posterior means of random effects}
#'    \item{model}{the fitted model}
#'    \item{acceptance.rates}{the achieved acceptance rates when using Metropolis-Hastings}
#'    \item{HPD.coverage}{coverage probability of HPD intervals (based on input)}
#'    \item{runtime}{total runtime of the sampler (measured in seconds)}
#'    \item{WAIC}{the Widely Applicable Information Criterion for model comparison
#'     (note that the WAIC is only computed by using the actually observed data,
#'     i.e., missing data are fully ignored when computing WAIC)}
#'    \item{learning.settings}{information on which parameters have been learned}
#'    \item{mcmc.settings}{details on the MCMC sampler}
#'  }
#' @export
#'
#' @references
#'
#'  Bitto, A. and Frühwirth-Schnatter, S. (2019). Achieving Shrinkage in a
#'  Time-Varying Parameter Model Framework. In: Journal of Econometrics, 210,
#'  75-97.
#'
#'  Frühwirth-Schnatter, S. and Wagner, H. (2010). Stochastic Model Specification
#'  Search for Gaussian and Partially Non-Gaussian State Space Models. Journal
#'  of Econometrics, 154, 85-100.
#'
#'  Pillow, J. and Scott, J. (2012). Fully Bayesian inference for neural models
#'  with negative-binomial spiking. In: Advances in neural information processing
#'  systems, 25.
#'
#'  Polson, N.G., Scott, J.G. and Windle, J. (2013). Bayesian Inference for
#'  Logistic Models using Pólya-Gamma Latent Variables. Journal of the American
#'  Statistical Association, 108, 1339-1349.
#'
#'  Pfeiler, R. and Wagner, H. (2024). Shrinkage in a Bayesian Panel Data Model
#'  with Time-Varying Coefficients. In: Einbeck, J. Maeng, H., Ogundium, E. and
#'  Perrakis, K. (Editors): Developments in Statistical Modelling, Springer, 109-115.
#'
#'
#'
#' @import stats
#'
#' @examples 1
panelTVP <- function(formula,
                     data,
                     model,
                     prior.reg = list(
                       d.tau = 0.001, e.tau = 0.001, d.xi = 0.001, e.xi = 0.001,
                       b.tau = 10, nu.tau = 5, b.xi = 10, nu.xi = 5,
                       a.tau = 1, kappa.tau = 10, a.xi = 1, kappa.xi = 10,
                       iota.tau = 1, iota.xi = 1,
                       learn.a.tau = TRUE, learn.a.xi = TRUE,
                       target.rate.tau = 0.44, target.rate.xi = 0.44,
                       learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                       type = "rw2", c = 1, B0 = 1
                     ),
                     prior.var = list(
                       learn.C0.hyp = list(g0 = 5, G0 = 3.333333), c0 = 2.5
                     ),
                     prior.load = list(
                       d.phi = 0.001, e.phi = 0.001, d.zeta = 0.001, e.zeta = 0.001,
                       a.phi = 1, kappa.phi = 10, a.zeta = 1, kappa.zeta = 10,
                       learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                       type = "rw2", c = 1, L0 = 1
                     ),
                     mcmc.opt = list(
                       chain.length = 12000, burnin = 2000, thin = 10, asis = TRUE
                     ),
                     settings.NegBin = list(
                       alpha.r = 2, beta.r = 1, expansion.steps = 10,
                       width = 1, p.overrelax = 0, accuracy.overrelax = 10
                     ),
                     HPD.coverage = 0.95
){

  # Initialization -------------------------------------------------------------

  if(prior.load$type == "cps"){
    tv.load = FALSE
  } else{
    tv.load = TRUE
  }

  miss <- ifelse(is.na(data$y), TRUE, FALSE)
  N.miss <- sum(miss)
  if(model == "Gaussian") data$y[miss] <- rnorm(n = N.miss)
  if(model %in% c("Probit", "Logit")) data$y[miss] <- rbinom(n = N.miss, size = 1, prob = 0.5)
  if(model == "NegBin") data$y[miss] <- MASS::rnegbin(n = N.miss, mu = 1, theta = 1)

  mf <- model.frame(formula = formula, data = data, drop.unused.levels = TRUE)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)

  tind <- data$t # here we need a check that t and id have to be columns in the data set !!!
  Tmax <- max(tind)
  id <- data$id

  y <- y[order(tind, id)]
  x <- x[order(tind, id), , drop = FALSE]
  df <- data.frame(tind, id)
  df <- df[order(tind, id),]

  df <- list(y = y, X = x, Tmax = Tmax,
             n = length(y)/max(tind),
             size = length(y), d = ncol(x),
             timeidx = df$tind, idx = df$id)

  alpha <- rnorm(df$d*2)
  sigma2 <- 1
  sigma2v <- sigma2

  C0 <- 1

  fi <- rep(0, df$n)
  fv <- rep(fi, df$Tmax)

  if(!tv.load){
    lambda <- 0
    reff <- lambda*fv
  } else{
    lambda <- rep(0, df$Tmax)
    reff <- c(t(matrix(lambda, ncol=df$n, nrow=df$Tmax)))*fv
  }

  alpha_lambda <- matrix(c(1.2,0.5))

  if(tv.load & length(prior.load$L0) == 1){
    prior.load$L0 <- rep(prior.load$L0, Tmax)
  }

  prior.reg$tau <- rep(10, df$d)
  prior.reg$xi <- rep(10, df$d)
  prior.load$phi <- 1
  prior.load$zeta <- 1

  # create return matrix for the MCMC samples

  if(prior.reg$type == "ind"){
    namesbetat <- unlist(lapply(1:df$d, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat <- unlist(lapply(1:df$d, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta <- paste0("beta", 1:df$d)
    namestheta <- paste0("theta", 1:df$d)
    namestau <- paste0("tau2", 1:df$d)
    namesxi <- paste0("xi2", 1:df$d)
  }
  namessgma2 <- "sigma2"
  if(!tv.load){
    nameslambdat <-"lambda_t"
  } else{
    nameslambdat <- paste0("lambda_t",1:df$Tmax)
  }
  if(prior.reg$type == "ind"){
    cnames <- c("SimNr", namesbetat, namessgma2, nameslambdat)
  } else{
    cnames <- c("SimNr",namesbetat, namesbeta, namestheta,
                namestau, namesxi, c("a.tau","kappa.tau","a.xi","kappa.xi"),
                namessgma2, nameslambdat)
  }
  if(prior.load$type=="rw1" | prior.load$type=="rw2"){
    cnames <- c(cnames,"lambda","psi","phi2","zeta2", "a.phi", "kappa.phi", "a.zeta", "kappa.zeta")
  }
  col_res <- length(cnames)
  res_frame <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res)
  colnames(res_frame) <- cnames

  f_sum <- rep(0, df$n)
  f_mat <- matrix(NA, nrow = (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin, ncol = length(y))

  # modification for Negative-Binomial model as we also want chain for r

  if(model == "NegBin"){
    cnames <- append(cnames, "r")
    col_res <- length(cnames)
    res_frame <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res)
    colnames(res_frame) <- cnames
  }

  # initialize acceptance rates

  ## r (negative binomial dispersion parameter)
  settings.NegBin$r.accept <- c()
  settings.NegBin$r.accept[1] <- 1 # we let metropolis start in 2nd iteration

  ## regression part
  prior.reg$xi.accept <- c()
  prior.reg$xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg$tau.accept <- c()
  prior.reg$tau.accept[1] <- 1 # we let metropolis start in 2nd iteration

  # fitting the model ----------------------------------------------------------

  if(model == "Gaussian"){
    result <- GaussianTVP(df = df,
                          prior.reg = prior.reg,
                          prior.var = prior.var,
                          prior.load = prior.load,
                          mcmc.opt = mcmc.opt,
                          sigma2v = sigma2v,
                          alpha = alpha,
                          lambda = lambda,
                          alpha_lambda = alpha_lambda,
                          reff = reff,
                          C0 = C0,
                          tv.load = tv.load,
                          res_frame = res_frame,
                          f_sum = f_sum,
                          f_mat = f_mat,
                          miss = miss,
                          HPD.coverage = HPD.coverage)
    class(result) <- "panelTVP.Gaussian"
  }

  if(model == "Probit"){
    result <- ProbitTVP(df = df,
                        prior.reg = prior.reg,
                        prior.load = prior.load,
                        mcmc.opt = mcmc.opt,
                        alpha = alpha,
                        lambda = lambda,
                        alpha_lambda = alpha_lambda,
                        reff = reff,
                        tv.load = tv.load,
                        res_frame = res_frame,
                        f_sum = f_sum,
                        f_mat = f_mat,
                        miss = miss,
                        HPD.coverage = HPD.coverage)
    class(result) <- "panelTVP.Probit"
  }

  if(model == "Logit"){
    result <- LogitTVP(df = df,
                       prior.reg = prior.reg,
                       prior.load = prior.load,
                       mcmc.opt = mcmc.opt,
                       alpha = alpha,
                       lambda = lambda,
                       alpha_lambda = alpha_lambda,
                       reff = reff,
                       tv.load = tv.load,
                       res_frame = res_frame,
                       f_sum = f_sum,
                       f_mat = f_mat,
                       miss = miss,
                       HPD.coverage = HPD.coverage)
    class(result) <- "panelTVP.Logit"
  }

  if(model == "NegBin"){
    result <- NegBinTVP(df = df,
                        prior.reg = prior.reg,
                        prior.load = prior.load,
                        mcmc.opt = mcmc.opt,
                        settings.NegBin = settings.NegBin,
                        alpha = alpha,
                        lambda = lambda,
                        alpha_lambda = alpha_lambda,
                        reff = reff,
                        tv.load = tv.load,
                        res_frame = res_frame,
                        f_sum = f_sum,
                        f_mat = f_mat,
                        HPD.coverage = HPD.coverage)
    class(result) <- "panelTVP.NegBin"
  }

  # add WAIC and remove chain of factor scores to save memory
  result$WAIC <- compute.waic(result)
  result[["fmcmc"]] <- NULL

  # adding learning settings to output
  hyperpara <- c("a.xi", "a.tau", "kappa.xi", "kappa.tau", "kappa.zeta", "kappa.phi")
  part <- c(rep("regression part", 4), rep("factor part", 2))
  learn <- c(prior.reg$learn.a.xi, prior.reg$learn.a.tau,
             prior.reg$learn.kappa.xi, prior.reg$learn.kappa.tau,
             prior.load$learn.kappa.zeta, prior.load$learn.kappa.phi)
  if(!(prior.reg$type %in% c("rw1", "rw2"))) learn[1:4] <- NA
  if(!(prior.load$type %in% c("rw1", "rw2"))) learn[5:6] <- NA
  result$learning.settings <- cbind(hyperpara, part, learn)
  colnames(result$learning.settings) <- c("hyperparameter", "model.part", "learned?")

  # adding mcmc setting to output (incl. ASIS Boolean)
  result$mcmc.settings <- mcmc.opt

  # rounding HPD lower bound to exactly 0 to cover cases that should be zero
  # but do not include zero due to sign flip
  index1 <- startsWith(rownames(result$posterior), "abs(")
  result$posterior[index1, "LO"] <- ifelse(result$posterior[index1,"LO"] < 0.01, 0,
                                           result$posterior[index1,"LO"])
  result$posterior["lambda_t1","LO"] <- ifelse(result$posterior["lambda_t1","LO"] < 0.01, 0,
                                               result$posterior["lambda_t1","LO"])

  return(result)

}
