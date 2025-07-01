#' Fit a Bayesian panel model with time-varying parameters for overdispersed and zero-inflated count data
#'
#' @param formula A formula argument consisting of two parts: the count model and
#'  the zero-inflation model (no default, see details)
#' @param data A data.frame containing the variables specified in the formula argument
#'  (no default)
#' @param prior.reg_nb A list of arguments for estimating the parameters of the regression
#'  part of the count model. The arguments are:
#'   \itemize{
#'    \item \code{d.tau}: shape parameter of Gamma prior for \eqn{\kappa^\tau} (count model)
#'    \item \code{e.tau}: rate parameter of Gamma prior for \eqn{\kappa^\tau} (count model)
#'    \item \code{d.xi}: shape parameter of Gamma prior for \eqn{\kappa^\xi} (count model)
#'    \item \code{e.xi}: rate parameter of Gamma prior for \eqn{\kappa^\xi} (count model)
#'    \item \code{b.tau}: part of the rate parameter of the Gamma prior for \eqn{a^\tau} (count model)
#'    \item \code{nu.tau}: shape parameter of the Gamma prior for \eqn{a^\tau} (count model)
#'    \item \code{b.xi}: part of the rate parameter of the Gamma prior for \eqn{a^\xi} (count model)
#'    \item \code{nu.xi}: shape parameter of the Gamma prior for \eqn{a^\xi} (count model)
#'    \item \code{a.tau}: shape parameter of the Gamma prior for \eqn{\tau^2_j} (count model)
#'    \item \code{kappa.tau}: part of the rate parameter of the Gamma prior for \eqn{\tau^2_j} (count model)
#'    \item \code{a.xi}: shape parameter of the Gamma prior for \eqn{\xi^2_j} (count model)
#'    \item \code{kappa.xi}: part of the rate parameter of the Gamma prior for \eqn{\xi^2_j} (count model)
#'    \item \code{iota.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\tau} (count model)
#'    \item \code{iota.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\xi} (count model)
#'    \item \code{learn.a.tau}: if TRUE \eqn{a^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.tau} used as starting value, if FALSE
#'     \eqn{a^\tau =} \code{a.tau} (count model)
#'    \item \code{learn.a.xi}: if TRUE \eqn{a^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.xi} used as starting value, if FALSE
#'     \eqn{a^\xi =} \code{a.xi} (count model)
#'    \item \code{tau.target.rate}: desired acceptance rate when updating
#'     \eqn{a^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.tau = FALSE}) (count model)
#'    \item \code{xi.target.rate}: desired acceptance rate when updating
#'     \eqn{a^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.xi = FALSE} (count model)
#'    \item \code{learn.kappa.tau}: if TRUE \eqn{\kappa^\tau} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.tau} used as starting value,
#'     if FALSE \eqn{\kappa^\tau = } \code{kappa.tau} (count model)
#'    \item \code{learn.kappa.xi}: if TRUE \eqn{\kappa^\xi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.xi} used as starting value,
#'     if FALSE \eqn{\kappa^\xi = } \code{kappa.xi} (count model)
#'    \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "rw1" (RW1 shrinkage prior), "rw2" (RW2 shrinkage prior)
#'      or "ind" (independence prior) (count model)
#'    \item \code{c}: prior parameter that scales the variance when using the
#'     random walk shrinkage prior (ignored when using independence prior) (count model)
#'    \item \code{B0}: prior variance on the regression parameters when using the
#'     independence prior (ignored when using shrinkage prior) (count model)
#'   }
#' @param prior.load_nb A list of arguments for estimating the parameters of the factor
#'  part of the count model. The arguments are:
#'  \itemize{
#'   \item \code{d.phi}: shape parameter of Gamma prior for \eqn{\kappa^\phi} (count model)
#'   \item \code{e.phi}: rate parameter of Gamma prior for \eqn{\kappa^\phi} (count model)
#'   \item \code{d.zeta}: shape parameter of Gamma prior for \eqn{\kappa^\zeta} (count model)
#'   \item \code{e.zeta}: rate parameter of Gamma prior for \eqn{\kappa^\zeta} (count model)
#'   \item \code{a.phi}: shape parameter of the Gamma prior for \eqn{\phi^2_j} (count model)
#'   \item \code{kappa.phi}: part of the rate parameter of the Gamma prior for \eqn{\phi^2_j} (count model)
#'   \item \code{a.zeta}: shape parameter of the Gamma prior for \eqn{\zeta^2_j} (count model)
#'   \item \code{kappa.zeta}: part of the rate parameter of the Gamma prior for \eqn{\zeta^2_j} (count model)
#'   \item \code{learn.kappa.phi}: if TRUE \eqn{\kappa^\phi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.phi} used as starting value,
#'     if FALSE \eqn{\kappa^\phi = } \code{kappa.phi} (count model)
#'   \item \code{learn.kappa_zeta}: if TRUE \eqn{\kappa^\zeta} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.zeta} used as starting value,
#'     if FALSE \eqn{\kappa^\zeta = } \code{kappa.zeta} (count model)
#'   \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "cps" (compound symmetric), "rw1" (RW1 shrinkage prior),
#'     "rw2" (RW2 shrinkage prior) or "ind" (independence prior) (count model)
#'   \item \code{c}: prior parameter that scales the variance when using the
#'     random walk shrinkage prior (ignored when using cps or independence prior) (count model)
#'   \item \code{L0}: prior variance on the factor loading (ignored when using shrinkage prior) (count model)
#'  }
#'  Note that for the factor model, the hyperparameters \code{a_phi} and
#'  \code{a_zeta} have to be fixed and are not sampled using Metropolis-Hastings
#' @param prior.reg_logit A list of arguments for estimating the parameters of the regression
#'  part of the zero-inflation model. The arguments are:
#'   \itemize{
#'    \item \code{d.tau}: shape parameter of Gamma prior for \eqn{\kappa^\tau} (zero-inflation model)
#'    \item \code{e.tau}: rate parameter of Gamma prior for \eqn{\kappa^\tau} (zero-inflation model)
#'    \item \code{d.xi}: shape parameter of Gamma prior for \eqn{\kappa^\xi} (zero-inflation model)
#'    \item \code{e.xi}: rate parameter of Gamma prior for \eqn{\kappa^\xi} (zero-inflation model)
#'    \item \code{b.tau}: part of the rate parameter of the Gamma prior for \eqn{a^\tau} (zero-inflation model)
#'    \item \code{nu.tau}: shape parameter of the Gamma prior for \eqn{a^\tau} (zero-inflation model)
#'    \item \code{b.xi}: part of the rate parameter of the Gamma prior for \eqn{a^\xi} (zero-inflation model)
#'    \item \code{nu.xi}: shape parameter of the Gamma prior for \eqn{a^\xi} (zero-inflation model)
#'    \item \code{a.tau}: shape parameter of the Gamma prior for \eqn{\tau^2_j} (zero-inflation model)
#'    \item \code{kappa.tau}: part of the rate parameter of the Gamma prior for \eqn{\tau^2_j} (zero-inflation model)
#'    \item \code{a.xi}: shape parameter of the Gamma prior for \eqn{\xi^2_j} (zero-inflation model)
#'    \item \code{kappa.xi}: part of the rate parameter of the Gamma prior for \eqn{\xi^2_j} (zero-inflation model)
#'    \item \code{iota.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\tau} (zero-inflation model)
#'    \item \code{iota.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\xi} (zero-inflation model)
#'    \item \code{learn.a.tau}: if TRUE \eqn{a^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.tau} used as starting value, if FALSE
#'     \eqn{a^\tau =} \code{a.tau} (zero-inflation model)
#'    \item \code{learn.a.xi}: if TRUE \eqn{a^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.xi} used as starting value, if FALSE
#'     \eqn{a^\xi =} \code{a.xi} (zero-inflation model)
#'    \item \code{tau.target.rate}: desired acceptance rate when updating
#'     \eqn{a^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.tau = FALSE}) (zero-inflation model)
#'    \item \code{xi.target.rate}: desired acceptance rate when updating
#'     \eqn{a^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.xi = FALSE} (zero-inflation model)
#'    \item \code{learn.kappa.tau}: if TRUE \eqn{\kappa^\tau} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.tau} used as starting value,
#'     if FALSE \eqn{\kappa^\tau = } \code{kappa.tau} (zero-inflation model)
#'    \item \code{learn.kappa.xi}: if TRUE \eqn{\kappa^\xi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.xi} used as starting value,
#'     if FALSE \eqn{\kappa^\xi = } \code{kappa.xi} (zero-inflation model)
#'    \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "rw1" (RW1 shrinkage prior), "rw2" (RW2 shrinkage prior)
#'      or "ind" (independence prior) (zero-inflation model)
#'    \item \code{c}: prior parameter that scales the variance when using the
#'     random walk shrinkage prior (ignored when using independence prior) (zero-inflation model)
#'    \item \code{B0}: prior variance on the regression parameters when using the
#'     independence prior (ignored when using shrinkage prior) (zero-inflation model)
#'   }
#' @param prior.load_logit A list of arguments for estimating the parameters of the factor
#'  part of the zero-inflation model. The arguments are:
#'  \itemize{
#'   \item \code{d.phi}: shape parameter of Gamma prior for \eqn{\kappa^\phi} (zero-inflation model)
#'   \item \code{e.phi}: rate parameter of Gamma prior for \eqn{\kappa^\phi} (zero-inflation model)
#'   \item \code{d.zeta}: shape parameter of Gamma prior for \eqn{\kappa^\zeta} (zero-inflation model)
#'   \item \code{e.zeta}: rate parameter of Gamma prior for \eqn{\kappa^\zeta} (zero-inflation model)
#'   \item \code{a.phi}: shape parameter of the Gamma prior for \eqn{\phi^2_j} (zero-inflation model)
#'   \item \code{kappa.phi}: part of the rate parameter of the Gamma prior for \eqn{\phi^2_j} (zero-inflation model)
#'   \item \code{a.zeta}: shape parameter of the Gamma prior for \eqn{\zeta^2_j} (zero-inflation model)
#'   \item \code{kappa.zeta}: part of the rate parameter of the Gamma prior for \eqn{\zeta^2_j} (zero-inflation model)
#'   \item \code{learn.kappa.phi}: if TRUE \eqn{\kappa^\phi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.phi} used as starting value,
#'     if FALSE \eqn{\kappa^\phi = } \code{kappa.phi} (zero-inflation model)
#'   \item \code{learn.kappa_zeta}: if TRUE \eqn{\kappa^\zeta} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.zeta} used as starting value,
#'     if FALSE \eqn{\kappa^\zeta = } \code{kappa.zeta} (zero-inflation model)
#'   \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "cps" (compound symmetric), "rw1" (RW1 shrinkage prior),
#'     "rw2" (RW2 shrinkage prior) or "ind" (independence prior) (zero-inflation model)
#'   \item \code{c}: prior parameter that scales the variance when using the
#'     random walk shrinkage prior (ignored when using cps or independence prior) (zero-inflation model)
#'   \item \code{L0}: prior variance on the factor loading (ignored when using shrinkage prior) (zero-inflation model)
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
#'  parameter r in count component of the ZINB model using univariate Slice sampling.
#'  The arguments are:
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
#' @param HPD.coverage Coverage probability of highest posterior density intervals
#'  (default yields 95 percent coverage)
#'
#' @description
#'   This is the main function for fitting a flexible Bayesian panel data model
#'   in the time-varying parameter framework for zero-inflated and overdispersed
#'   count outcomes. By using shrinkage priors, it is
#'   possible to identify whether an effect is time-varying, time-invariant or zero.
#'   If you want to fit a model to a Gaussian, binary or Negative Binomial
#'   response, you should use the function [panelTVP()] instead.
#'
#' @details
#'  This function fits a Bayesian time-varying parameter panel data model to a
#'  longitudinal response \eqn{y_{it}} for \eqn{i \in \{1,\dots,n\}} subjects that are observed
#'  at \eqn{t \in \{1,\dots,\text{T}\}}
#'  time point. By using this function, it is assumed that the outcome is a possibly overdispersed count
#'  variable subject to zero-inflation (i.e., more zeros are present than what is
#'  assumed by a standard count distribution). To account for zero-inflation,
#'  it is assumed that there are two latent classes of zeros: structural and
#'  at-risk zeros. A structural zero belongs to an observation that is not at risk of
#'  experiencing the event, whereas an at-risk zero belongs to an observation that is
#'  at-risk of experiencing the event but has for some reason a zero recorded.
#'  It is therefore assumed that the outcome \eqn{y_{it}} is a realization of a
#'  mixture distribution with a point mass at zero (for the structural zeros) and
#'  a standard Negative Binomial count model for observations that are at-risk,
#'  i.e., this includes both at-risk zeros and positive counts. The ZINB model
#'  can, thus, be stated as (Neelon, 2019)
#'  \deqn{ y_{it} \sim (1-\pi_{it}) \cdot \mathbb{I}_{(w_{it} = 0)} + \pi_{it} \cdot \mathcal{NB}(q_{it},r) \mathbb{I}_{(w_{it}=1)},}
#'  where \eqn{w_{it}} is a latent at-risk indicator such that with probability
#'  \eqn{1-\pi_{it}, w_{it} = 0} and with probability \eqn{\pi_{it}, w_{it} = 1}.
#'  The Negative Binomial component is parameterized such that its mass function
#'  is given by
#'  \deqn{p(y_{it}|r,q) = \frac{\Gamma(y_{it}+r)}{\Gamma(r)y_{it}!} (1-q_{it})^r q_{it}^{y_{it}},}
#'  where \eqn{r} is a common dispersion parameter. We use this parameterization
#'  of the Negative Binomial distribution
#'  to fascilitate posterior inference with Pólya-Gamma random variables
#'  (see Pillow and Scott, 2012).
#'  Moreover, there are two separate linear predictors in the ZINB model. The first
#'  linear predictor \eqn{\eta_{it}^{\text{logit}}} is related to the zero-inflated
#'  part of the model through
#'  \deqn{\pi_{it} = \frac{\exp(\eta_{it}^{\text{logit}})}{1+\exp(\eta_{it}^{\text{logit}})},}
#'  whereas the second linear predictor \eqn{\eta_{it}^{\text{nb}}} is related to
#'  the count part of the model though
#'  \deqn{q_{it} = \frac{\exp(\eta_{it}^{\text{nb}})}{1+\exp(\eta_{it}^{\text{nb}})}.}
#'  In every iteration of the MCMC sampler, the latent at-risk indicators \eqn{w_{it}} are updated
#'  and the parameters in the count model are estimated by using only the observations
#'  that are currently in the risk set. Note that for \eqn{y_{it} > 0} the at-risk indicators
#'  are fixed at \eqn{w_{it} = 1} in every iteration.

#'  In the following, we describe the general
#'  idea of the model. To simplify notation, we suppress the
#'  superscript and use \eqn{\eta_{it}} for both \eqn{\eta_{it}^{\text{logit}}} and
#'  \eqn{\eta_{it}^{\text{nb}}}. However, note that different sets of covariates
#'  can enter the two predictors.
#'
#'  The function \code{panelTVP_ZINB} allows for time-varying parameters in both
#'  linear predictors. To avoid overfitting, we place shrinkage priors on all parameters
#'  of the model. To facilitate this, we consider the non-centered parameterization
#'  of the model (see Frühwirth-Schnatter and Wagner (2010) for details on this
#'  parameterization), where the linear predictors are given by
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
#' @returns The function returns an object of class \code{panelTVP.ZINB}.
#'  The returned object contains a list of the following elements
#'  \describe{
#'    \item{data}{the data used for fitting the model and additional context information
#'    derived from the data}
#'    \item{mcmc_logit}{Markov Chains for every parameter except for the factor scores
#'     (to save memory) for the Logit component of the model}
#'    \item{mcmc_nb}{Markov Chains for every parameter except for the factor scores
#'     (to save memory) for the Negative Binomial component of the model}
#'    \item{posterior_logit}{preliminary summary of posterior results for the
#'     Logit component of the model}
#'    \item{posterior_nb}{preliminary summary of posterior results for the
#'     Negative Binomial component of the model}
#'    \item{fmean_logit}{posterior means of random effects for the
#'      Logit component of the model}
#'    \item{fmean_nb}{posterior means of random effects for the
#'      Negative Binomial component of the model}
#'    \item{model}{the fitted model (here only Bayesian Zero-Inflated Negative Binomial possible)}
#'    \item{acceptance.rates}{the achieved acceptance rates when using
#'       Metropolis-Hastings for both components of the model}
#'    \item{HPD.coverage}{coverage probability of HPD intervals (based on input)}
#'    \item{runtime}{total runtime of the sampler (measured in seconds)}
#'    \item{WAIC}{the Widely Applicable Information Criterion for model comparison}
#'    \item{learning.settings_logit}{information on which parameters have been learned
#'       for the Logit component of the model}
#'    \item{learning.settings_nb}{information on which parameters have been learned
#'       for the Negative Binomial component of the model}
#'    \item{mcmc.settings}{details on general MCMC settings}
#'  }
#'
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
#'  Neelon, B. (2019). Bayesian Zero-Inflated Negative Binomial Regression Based
#'  on Pólya-Gamma Mixtures. In: Bayesian Analysis, 14, 829-855.
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
#'
#' @examples 0
panelTVP_ZINB <- function(formula,
                          data,
                          prior.reg_nb = list(
                            d.tau = 0.001, e.tau = 0.001, d.xi = 0.001, e.xi = 0.001,
                            b.tau = 10, nu.tau = 5, b.xi = 10, nu.xi = 5,
                            a.tau = 1, kappa.tau = 10, a.xi = 1, kappa.xi = 10,
                            iota.tau = 1, iota.xi = 1,
                            learn.a.tau = TRUE, learn.a.xi = TRUE,
                            target.rate.tau = 0.44, target.rate.xi = 0.44,
                            learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                            type = "rw2", c = 1, B0 = 1
                          ),
                          prior.load_nb = list(
                            d.phi = 0.001, e.phi = 0.001, d.zeta = 0.001, e.zeta = 0.001,
                            a.phi = 1, kappa.phi = 10, a.zeta = 1, kappa.zeta = 10,
                            learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                            type = "rw2", c = 1, L0 = 1
                          ),
                          prior.reg_logit = list(
                            d.tau = 0.001, e.tau = 0.001, d.xi = 0.001, e.xi = 0.001,
                            b.tau = 10, nu.tau = 5, b.xi = 10, nu.xi = 5,
                            a.tau = 1, kappa.tau = 10, a.xi = 1, kappa.xi = 10,
                            iota.tau = 1, iota.xi = 1,
                            learn.a.tau = TRUE, learn.a.xi = TRUE,
                            target.rate.tau = 0.44, target.rate.xi = 0.44,
                            learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                            type = "rw2", c = 1, B0 = 1
                          ),
                          prior.load_logit = list(
                            d.phi = 0.001, e.phi = 0.001, d.zeta = 0.001, e.zeta = 0.001,
                            a.phi = 1, kappa.phi = 10, a.zeta = 1, kappa.zeta = 10,
                            learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                            type = "rw2", c = 1, L0 = 1
                          ),
                          mcmc.opt = list(
                            chain.length = 12000, thin = 10, burnin = 2000, asis = TRUE
                          ),
                          settings.NegBin = list(
                            alpha.r = 2, beta.r = 1, expansion.steps = 20,
                            width = 1, p.overrelax = 0, accuracy.overrelax = 20
                          ),
                          HPD.coverage = 0.95
){

  # Initialization

  if(prior.load_nb$type == "cps"){
    tv.load_nb = FALSE
  } else{
    tv.load_nb = TRUE
  }

  if(prior.load_logit$type == "cps"){
    tv.load_logit = FALSE
  } else{
    tv.load_logit = TRUE
  }

  miss <- ifelse(is.na(data$y), TRUE, FALSE)
  N.miss <- sum(miss)
  data$y[miss] <- MASS::rnegbin(N.miss, mu = 1, theta = 2)

  rhs <- formula[[3]]
  rhs_str <- paste(deparse(rhs), collapse = "")
  rhs.parts <- strsplit(rhs_str, "\\|")[[1]]
  rhs.parts <- trimws(rhs.parts)
  formula_nb <- as.formula(paste(deparse(formula[[2]]), "~", rhs.parts[1]))
  formula_logit  <- as.formula(paste(deparse(formula[[2]]), "~", rhs.parts[2]))
  mf_nb <- model.frame(formula = formula_nb, data = data, drop.unused.levels = TRUE)
  mt_nb <- attr(mf_nb, "terms")
  x_nb <- model.matrix(mt_nb, mf_nb)
  mf_logit <- model.frame(formula = formula_logit, data = data, drop.unused.levels = TRUE)
  mt_logit <- attr(mf_logit, "terms")
  x_logit <- model.matrix(mt_logit, mf_logit)
  y <- model.response(mf_nb)
  tind <- data$t
  Tmax <- max(tind)
  id <- data$id
  y <- y[order(tind, id)]
  x_nb <- x_nb[order(tind, id), , drop = FALSE]
  x_logit <- x_logit[order(tind, id), , drop = FALSE]
  df <- data.frame(tind, id)
  df <- df[order(tind, id),]
  df <- list(y = y, X_nb = x_nb, X_logit = x_logit,
             Tmax = Tmax, n = length(y)/max(tind), size = length(y),
             d_nb = ncol(x_nb), d_logit = ncol(x_logit),
             timeidx = df$tind, idx = df$id)

  alpha_nb <- rnorm(df$d_nb*2)
  alpha_logit <- rnorm(df$d_logit*2)

  fi_nb <- rep(0, df$n)
  fv_nb <- rep(fi_nb, df$Tmax)
  fi_logit <- rep(0, df$n)
  fv_logit <- rep(fi_logit, df$Tmax)

  if(!tv.load_nb){
    lambda_nb <- 0
    reff_nb <- lambda_nb*fv_nb
  } else{
    lambda_nb <- rep(0, df$Tmax)
    reff_nb <- c(t(matrix(lambda_nb, ncol=df$n, nrow=df$Tmax)))*fv_nb
  }
  if(!tv.load_logit){
    lambda_logit <- 0
    reff_logit <- lambda_logit*fv_logit
  } else{
    lambda_logit <- rep(0, df$Tmax)
    reff_logit <- c(t(matrix(lambda_logit, ncol=df$n, nrow=df$Tmax)))*fv_logit
  }

  alpha_lambda_nb <- matrix(c(1.2,0.5))
  alpha_lambda_logit <- matrix(c(1.2,0.5))

  if(tv.load_nb & length(prior.load_nb$L0) == 1){
    prior.load_nb$L0 <- rep(prior.load_nb$L0, Tmax)
  }
  if(tv.load_logit & length(prior.load_logit$L0) == 1){
    prior.load_logit$L0 <- rep(prior.load_logit$L0, Tmax)
  }

  prior.reg_nb$tau <- rep(10, df$d_nb)
  prior.reg_nb$xi <- rep(10, df$d_nb)
  prior.load_nb$phi <- 1
  prior.load_nb$zeta <- 1

  prior.reg_logit$tau <- rep(10, df$d_logit)
  prior.reg_logit$xi <- rep(10, df$d_logit)
  prior.load_logit$phi <- 1
  prior.load_logit$zeta <- 1

  #-----------------------------------------------------------------------------
  # create return matrix for the MCMC samples

  if(prior.reg_nb$type == "ind"){
    namesbetat_nb <- unlist(lapply(1:df$d_nb, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat_nb <- unlist(lapply(1:df$d_nb, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta_nb <- paste0("beta", 1:df$d_nb)
    namestheta_nb <- paste0("theta", 1:df$d_nb)
    namestau_nb <- paste0("tau2", 1:df$d_nb)
    namesxi_nb <- paste0("xi2", 1:df$d_nb)
  }
  if(prior.reg_logit$type == "ind"){
    namesbetat_logit <- unlist(lapply(1:df$d_logit, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat_logit <- unlist(lapply(1:df$d_logit, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta_logit <- paste0("beta", 1:df$d_logit)
    namestheta_logit <- paste0("theta", 1:df$d_logit)
    namestau_logit <- paste0("tau2", 1:df$d_logit)
    namesxi_logit <- paste0("xi2", 1:df$d_logit)
  }

  if(!tv.load_nb){
    nameslambdat_nb <-"lambda_t"
  } else{
    nameslambdat_nb <- paste0("lambda_t",1:df$Tmax)
  }
  if(prior.reg_nb$type == "ind"){
    cnames_nb <- c("SimNr", namesbetat_nb, nameslambdat_nb)
  } else{
    cnames_nb <- c("SimNr",namesbetat_nb, namesbeta_nb, namestheta_nb,
                   namestau_nb, namesxi_nb, c("a.tau","kappa.tau","a.xi","kappa.xi"),
                   nameslambdat_nb)
  }
  if(prior.load_nb$type=="rw1" | prior.load_nb$type=="rw2"){
    cnames_nb <- c(cnames_nb,"lambda","psi","phi2","zeta2", "a.phi", "kappa.phi",
                   "a.zeta", "kappa.zeta")
  }
  cnames_nb <- append(cnames_nb, "r")
  col_res_nb <- length(cnames_nb)
  res_frame_nb <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res_nb)
  colnames(res_frame_nb) <- cnames_nb
  if(!tv.load_logit){
    nameslambdat_logit <-"lambda_t"
  } else{
    nameslambdat_logit <- paste0("lambda_t",1:df$Tmax)
  }
  if(prior.reg_logit$type == "ind"){
    cnames_logit <- c("SimNr", namesbetat_logit, nameslambdat_logit)
  } else{
    cnames_logit <- c("SimNr",namesbetat_logit, namesbeta_logit, namestheta_logit,
                      namestau_logit, namesxi_logit, c("a.tau","kappa.tau","a.xi","kappa.xi"),
                      nameslambdat_logit)
  }
  col_res_logit <- length(cnames_logit)
  if(prior.load_logit$type=="rw1" | prior.load_logit$type=="rw2"){
    cnames_logit <- c(cnames_logit,"lambda","psi","phi2","zeta2", "a.phi", "kappa.phi",
                      "a.zeta", "kappa.zeta")
    col_res_logit <- length(cnames_logit)
  }
  res_frame_logit <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res_logit)
  colnames(res_frame_logit) <- cnames_logit

  f_sum_nb <- rep(0, df$n)
  f_mat_nb <- matrix(NA, nrow = (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin,
                     ncol = length(y))
  f_sum_logit <- rep(0, df$n)
  f_mat_logit <- matrix(NA, nrow = (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin,
                        ncol = length(y))

  # initialize acceptance rate for dispersion parameter r
  settings.NegBin$r.accept <- c()

  ## regression part
  prior.reg_nb$xi.accept <- c()
  prior.reg_nb$xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg_nb$tau.accept <- c()
  prior.reg_nb$tau.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg_logit$xi.accept <- c()
  prior.reg_logit$xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg_logit$tau.accept <- c()
  prior.reg_logit$tau.accept[1] <- 1 # we let metropolis start in 2nd iteration

  # Fitting the ZINB model

  result <- zinbTVP(df = df,
                    prior.reg_nb = prior.reg_nb,
                    prior.reg_logit = prior.reg_logit,
                    prior.load_nb = prior.load_nb,
                    prior.load_logit = prior.load_logit,
                    mcmc.opt = mcmc.opt,
                    settings.NegBin = settings.NegBin,
                    alpha_nb = alpha_nb,
                    alpha_logit = alpha_logit,
                    lambda_nb = lambda_nb,
                    lambda_logit = lambda_logit,
                    alpha_lambda_nb = alpha_lambda_nb,
                    alpha_lambda_logit = alpha_lambda_logit,
                    reff_nb = reff_nb,
                    reff_logit = reff_logit,
                    tv.load_nb = tv.load_nb,
                    tv.load_logit = tv.load_logit,
                    res_frame_nb = res_frame_nb,
                    res_frame_logit = res_frame_logit,
                    f_sum_nb = f_sum_nb,
                    f_sum_logit = f_sum_logit,
                    f_mat_nb = f_mat_nb,
                    f_mat_logit = f_mat_logit,
                    miss = miss,
                    HPD.coverage = HPD.coverage)
  class(result) <- "panelTVP.ZINB"

  # add WAIC and remove chain of factor scores and risk-indicators to save memory
  result$WAIC <- compute.waic(result)
  result[["fmcmc_logit"]] <- NULL
  result[["fmcmc_nb"]] <- NULL
  result[["mcmc_risk"]] <- NULL

  # adding learning settings to output
  hyperpara <- c("a.xi", "a.tau", "kappa.xi", "kappa.tau", "kappa.zeta", "kappa.phi")
  part <- c(rep("regression part", 4), rep("factor part", 2))
  learn_nb <- c(prior.reg_nb$learn.a.xi, prior.reg_nb$learn.a.tau,
                prior.reg_nb$learn.kappa.xi, prior.reg_nb$learn.kappa.tau,
                prior.load_nb$learn.kappa.zeta, prior.load_nb$learn.kappa.phi)
  if(!(prior.reg_nb$type %in% c("rw1", "rw2"))) learn_nb[1:4] <- NA
  if(!(prior.load_nb$type %in% c("rw1", "rw2"))) learn_nb[5:6] <- NA
  learn_logit <- c(prior.reg_logit$learn.a.xi, prior.reg_logit$learn.a.tau,
                   prior.reg_logit$learn.kappa.xi, prior.reg_logit$learn.kappa.tau,
                   prior.load_logit$learn.kappa.zeta, prior.load_logit$learn.kappa.phi)
  if(!(prior.reg_logit$type %in% c("rw1", "rw2"))) learn_logit[1:4] <- NA
  if(!(prior.load_logit$type %in% c("rw1", "rw2"))) learn_logit[5:6] <- NA
  result$learning.settings_nb <- cbind(hyperpara, part, learn_nb)
  colnames(result$learning.settings_nb) <- c("hyperparameter", "model.part", "learned?")
  result$learning.settings_logit <- cbind(hyperpara, part, learn_logit)
  colnames(result$learning.settings_logit) <- c("hyperparameter", "model.part", "learned?")

  # adding mcmc setting to output (incl. ASIS Boolean)
  result$mcmc.settings <- mcmc.opt

  # rounding HPD lower bound to exactly 0 to cover cases that should be zero
  # but do not include zero due to sign flip
  index1_nb <- startsWith(rownames(result$posterior_nb), "abs(")
  result$posterior_nb[index1_nb, "LO"] <- ifelse(result$posterior_nb[index1_nb,"LO"] < 0.01, 0,
                                                 result$posterior_nb[index1_nb,"LO"])
  result$posterior_nb["lambda_t1","LO"] <- ifelse(result$posterior_nb["lambda_t1","LO"] < 0.01, 0,
                                                  result$posterior_nb["lambda_t1","LO"])
  index1_logit <- startsWith(rownames(result$posterior_logit), "abs(")
  result$posterior_logit[index1_logit, "LO"] <- ifelse(result$posterior_logit[index1_logit,"LO"] < 0.01, 0,
                                                       result$posterior_logit[index1_logit,"LO"])
  result$posterior_logit["lambda_t1","LO"] <- ifelse(result$posterior_logit["lambda_t1","LO"] < 0.01, 0,
                                                     result$posterior_logit["lambda_t1","LO"])

  return(result)

}

