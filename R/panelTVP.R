#' Fit a Bayesian panel model with time-varying parameters
#'
#' @usage panelTVP(formula = NULL,
#'          data = NULL,
#'          id = NULL,
#'          t = NULL,
#'          model = NULL,
#'          prior.reg = list(),
#'          prior.var = list(),
#'          prior.load = list(),
#'          prior.reg_zinb.count = list(),
#'          prior.load_zinb.count = list(),
#'          prior.reg_zinb.inflation = list(),
#'          prior.load_zinb.inflation = list(),
#'          mcmc.opt = list(),
#'          settings.NegBin = list(),
#'          HPD.coverage = 0.95,
#'          R.WAIC = 20,
#'          posterior.predictive.matrix = FALSE,
#'          random.effects = TRUE,
#'          progress.bar = FALSE)
#'
#' @param formula the usual formula argument in regression methods, e.g., as in [lm()].
#'  When fitting a Zero-Inflated Negative Binomial model, the covariates for the
#'  Negative Binomial (count) component and the Logit (zero-inflation) component are
#'  separated, e.g., \code{y ~ W1.nb | W1.logit + W2.logit} will consider the variable \code{W1.nb}
#'  for the count component and the variables \code{W1.logit} as well as \code{W2.logit} for the
#'  zero-inflation component of the model (no default)
#' @param data a data frame that contains the variables of the formula argument
#'  (no default)
#' @param id a vector with length equal to the number of observations in the dataset, i.e.,
#'  this is the 'subject' variable in the dataset (no default)
#' @param t a vector with length equal to the number of observations in the dataset, i.e.,
#'  this is the 'time' variable in the dataset (no default)
#' @param model a character indicating which model should be estimated.
#'  This parameter is either 'Gaussian', 'Probit', 'Logit', 'NegBin' or 'ZINB' depending on
#'  whether a model for Gaussian, Probit, Logit, Negative Binomial or
#'  Zero-Inflated Negative Binomial response data should be fitted (no default)
#' @param prior.reg a list of arguments for estimating the parameters of the regression
#'  part of the model. This argument is ignored when \code{model = 'ZINB'}. The arguments are:
#'   \itemize{
#'    \item \code{a.tau}: hyperparameter for learning \eqn{\tau^2_j}
#'     (double Gamma, triple Gamma)
#'    \item \code{a.xi}: hyperparameter for learning \eqn{\xi^2_j}
#'     (double Gamma, triple Gamma)
#'    \item \code{learn.a.tau}: if TRUE \eqn{a^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.tau} used as starting value, if FALSE
#'     \eqn{a^\tau =} \code{a.tau} (not learned under alternative triple Gamma)
#'    \item \code{learn.a.xi}: if TRUE \eqn{a^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.xi} used as starting value, if FALSE
#'     \eqn{a^\xi =} \code{a.xi} (not learned under alternative triple Gamma)
#'     \item \code{alpha.a.tau}: hyperparameter for learning \eqn{a^\tau}
#'      (argument is ignored when \code{learn.a.tau = FALSE})
#'     \item \code{alpha.a.xi}: hyperparameter for learning \eqn{a^\xi}
#'      (argument is ignored when \code{learn.a.xi = FALSE})
#'     \item \code{beta.a.tau}: hyperparameter for learning \eqn{a^\tau}
#'      (argument is ignored when \code{learn.a.tau = FALSE})
#'     \item \code{beta.a.xi}: hyperparameter for learning \eqn{a^\xi}
#'      (argument is ignored when \code{learn.a.xi = FALSE})
#'    \item \code{iota.a.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\tau} (argument is ignored when
#'     \code{learn.a.tau = FALSE})
#'    \item \code{iota.a.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\xi} (argument is ignored when
#'     \code{learn.a.xi = FALSE})
#'    \item \code{target.rate.a.tau}: desired acceptance rate when updating
#'     \eqn{a^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.tau = FALSE})
#'    \item \code{target.rate.a.xi}: desired acceptance rate when updating
#'     \eqn{a^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.xi = FALSE})
#'    \item \code{c.tau:} shape parameter in the Gamma prior for \eqn{\kappa^\tau_j}
#'     (only triple Gamma)
#'    \item \code{c.xi:} shape parameter in the Gamma prior for \eqn{\kappa^\xi_j}
#'     (only triple Gamma)
#'    \item \code{learn.c.tau}: if TRUE \eqn{c^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{c.tau} used as starting value, if FALSE
#'     \eqn{c^\tau =} \code{c.tau} (not learned under alternative triple Gamma)
#'    \item \code{learn.c.xi}: if TRUE \eqn{c^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{c.xi} used as starting value, if FALSE
#'     \eqn{c^\xi =} \code{c.xi} (not learned under alternative triple Gamma)
#'     \item \code{alpha.c.tau}: hyperparameter for learning \eqn{c^\tau}
#'      (argument is ignored when \code{learn.c.tau = FALSE})
#'     \item \code{alpha.c.xi}: hyperparameter for learning \eqn{c^\xi}
#'      (argument is ignored when \code{learn.c.xi = FALSE})
#'     \item \code{beta.c.tau}: hyperparameter for learning \eqn{c^\tau}
#'      (argument is ignored when \code{learn.c.tau = FALSE})
#'     \item \code{beta.c.xi}: hyperparameter for learning \eqn{c^\xi}
#'      (argument is ignored when \code{learn.c.xi = FALSE})
#'    \item \code{iota.c.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{c^\tau}
#'    \item \code{iota.c.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{c^\xi}
#'    \item \code{target.rate.c.tau}: desired acceptance rate when updating
#'     \eqn{c^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.c.tau = FALSE})
#'    \item \code{target.rate.c.xi}: desired acceptance rate when updating
#'     \eqn{c^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.c.xi = FALSE})
#'    \item \code{kappa.tau}: hyperparameter for learning \eqn{\tau^2_j}
#'     (double Gamma, original triple Gamma). Note that this is the global parameter.
#'     The component-specific parameters are learned automatically when triple
#'     Gamma prior is selected
#'    \item \code{kappa.xi}: hyperparameter for learning \eqn{\xi^2_j}
#'     (double Gamma, original triple Gamma). Note that this is the global parameter.
#'     The component-specific parameters are learned automatically when triple
#'     Gamma prior is selected
#'    \item \code{learn.kappa.tau}: if TRUE \eqn{\kappa^\tau} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.tau} used as starting value,
#'     if FALSE \eqn{\kappa^\tau = } \code{kappa.tau}
#'    \item \code{learn.kappa.xi}: if TRUE \eqn{\kappa^\xi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.xi} used as starting value,
#'     if FALSE \eqn{\kappa^\xi = } \code{kappa.xi}
#'    \item \code{d.tau}: hyperparameter for learning \eqn{\kappa^\tau}
#'     (only double Gamma)
#'    \item \code{e.tau}: hyperparameter for learning \eqn{\kappa^\tau}
#'     (double Gamma, original triple Gamma prior). Note that under triple Gamma shrinkage
#'     this parameter is only the starting value as it is learned during MCMC
#'    \item \code{d.xi}: hyperparameter for learning \eqn{\kappa^\xi}
#'     (only double Gamma)
#'    \item \code{e.xi}: hyperparameter for learning \eqn{\kappa^\xi}
#'     (double Gamma, original triple Gamma). Note that under triple Gamma shrinkage
#'     this parameter is only the starting value as it is learned during MCMC
#'    \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "rw-t0" (shrinkage prior starting at T = 0),
#'     "rw-t1" (shrinkage prior starting at T = 1) or "ind" (independence prior)
#'    \item \code{c}: prior parameter that scales the variance when using
#'     shrinkage prior (ignored when using independence prior)
#'    \item \code{B0}: prior variance on the regression parameters when using
#'     independence prior (ignored when using shrinkage prior)
#'    \item \code{TG}: Boolean indicating whether triple Gamma prior should be used;
#'     when \code{TG = FALSE} and \code{type = "rw-t0/1"} then the double Gamma
#'     prior is used; when \code{type = "ind"} the argument \code{TG} is ignored
#'    \item \code{TG.alternative}: Boolean indicating whether the alternative
#'     triple Gamma representation of Knaus and Frühwirth-Schnatter (2025) should
#'     be used. If \code{TG.alternative = FALSE}, the original triple Gamma
#'     representation based on Cadonna et al. (2020) will be used instead.
#'     Note that the hyperparameters \eqn{a} and \eqn{c} will not be sampled
#'     and have to be set to fixed values in the alternative representation.
#'     If \code{TG = FALSE}, this argument is ignored
#'   }
#' @param prior.var a list of arguments for estimating the homoscedastic error variance
#'  in a Gaussian / Normal model. For other models, this argument is ignored.
#'  The arguments are:
#'  \itemize{
#'   \item \code{learn.C0.hyp}: this argument is a list containing the prior
#'    parameters for the prior rate \eqn{C_0} of \eqn{\sigma^2}. The parameters are:
#'    \itemize{
#'     \item \code{g0}: shape parameter of Gamma prior on \eqn{C_0}
#'     \item \code{G0}: rate parameter of Gamma prior on \eqn{C_0}
#'    }
#'   \item \code{c0}: shape parameter of Inverse-Gamma prior on \eqn{\sigma^2}
#'    }
#' @param prior.load a list of arguments for estimating the parameters of the factor
#'  part of the model.  This argument is ignored when \code{model = 'ZINB'}. The arguments are:
#'  \itemize{
#'   \item \code{a.phi}: hyperparameter for learning \eqn{\phi^2}
#'   \item \code{a.zeta}: hyperparameter for learning  \eqn{\zeta^2}
#'   \item \code{kappa.phi}: hyperparameter for learning \eqn{\phi^2}
#'   \item \code{kappa.zeta}: part of the rate parameter of the Gamma prior for \eqn{\zeta^2}
#'   \item \code{learn.kappa.phi}: if TRUE \eqn{\kappa^\phi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.phi} used as starting value,
#'     if FALSE \eqn{\kappa^\phi = } \code{kappa.phi}
#'   \item \code{learn.kappa.zeta}: if TRUE \eqn{\kappa^\zeta} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.zeta} used as starting value,
#'     if FALSE \eqn{\kappa^\zeta = } \code{kappa.zeta}
#'   \item \code{d.phi}: hyperparameter for learning \eqn{\kappa^\phi}
#'   \item \code{d.zeta}: hyperparameter for learning \eqn{\kappa^\zeta}
#'   \item \code{e.phi}: hyperparameter for learning \eqn{\kappa^\phi}
#'   \item \code{e.zeta}: hyperparameter for learning \eqn{\kappa^\zeta}
#'   \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "cps" (compound symmetric), "rw-t0" (shrinkage prior starting at T = 0),
#'     "rw-t1" (shrinkage prior starting at T = 1) or "ind" (independence prior)
#'   \item \code{c}: prior parameter that scales the variance when using shrinkage prior
#'    (ignored when using cps or independence prior)
#'   \item \code{L0}: prior variance on the factor loading when using either cps or independence
#'    prior (ignored when using shrinkage prior)
#'  }
#'  Note that for the factor model, the hyperparameters \code{a.phi} and
#'  \code{a.zeta} have to be fixed and are not sampled using Metropolis-Hastings.
#'  There is also no triple Gamma prior as only one factor is present
#' @param prior.reg_zinb.count A list of arguments for estimating the parameters of the regression
#'  part of the count model. This argument is only used when \code{model = 'ZINB'} and
#'  otherwise ignored. The arguments are:
#'   \itemize{
#'    \item \code{a.tau}: hyperparameter for learning \eqn{\tau^2_j}
#'     (double Gamma, triple Gamma) (count component)
#'    \item \code{a.xi}: hyperparameter for learning \eqn{\xi^2_j}
#'     (double Gamma, triple Gamma) (count component)
#'    \item \code{learn.a.tau}: if TRUE \eqn{a^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.tau} used as starting value, if FALSE
#'     \eqn{a^\tau =} \code{a.tau} (not learned under alternative triple Gamma)
#'     (count component)
#'    \item \code{learn.a.xi}: if TRUE \eqn{a^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.xi} used as starting value, if FALSE
#'     \eqn{a^\xi =} \code{a.xi} (not learned under alternative triple Gamma)
#'      (count component)
#'     \item \code{alpha.a.tau}: hyperparameter for learning \eqn{a^\tau}
#'      (argument is ignored when \code{learn.a.tau = FALSE}) (count component)
#'     \item \code{alpha.a.xi}: hyperparameter for learning \eqn{a^\xi}
#'      (argument is ignored when \code{learn.a.xi = FALSE}) (count component)
#'     \item \code{beta.a.tau}: hyperparameter for learning \eqn{a^\tau}
#'      (argument is ignored when \code{learn.a.tau = FALSE}) (count component)
#'     \item \code{beta.a.xi}: hyperparameter for learning \eqn{a^\xi}
#'      (argument is ignored when \code{learn.a.xi = FALSE}) (count component)
#'    \item \code{iota.a.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\tau} (argument is ignored when
#'     \code{learn.a.tau = FALSE}) (count component)
#'    \item \code{iota.a.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\xi} (argument is ignored when
#'     \code{learn.a.xi = FALSE}) (count component)
#'    \item \code{target.rate.a.tau}: desired acceptance rate when updating
#'     \eqn{a^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.tau = FALSE}) (count component)
#'    \item \code{target.rate.a.xi}: desired acceptance rate when updating
#'     \eqn{a^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.xi = FALSE}) (count component)
#'    \item \code{c.tau:} shape parameter in the Gamma prior for \eqn{\kappa^\tau_j}
#'     (only triple Gamma) (count component)
#'    \item \code{c.xi:} shape parameter in the Gamma prior for \eqn{\kappa^\xi_j}
#'     (only triple Gamma) (count component)
#'    \item \code{learn.c.tau}: if TRUE \eqn{c^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{c.tau} used as starting value, if FALSE
#'     \eqn{c^\tau =} \code{c.tau} (not learned under alternative triple Gamma)
#'      (count component)
#'    \item \code{learn.c.xi}: if TRUE \eqn{c^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{c.xi} used as starting value, if FALSE
#'     \eqn{c^\xi =} \code{c.xi} (not learned under alternative triple Gamma)
#'      (count component)
#'     \item \code{alpha.c.tau}: hyperparameter for learning \eqn{c^\tau}
#'      (argument is ignored when \code{learn.c.tau = FALSE}) (count component)
#'     \item \code{alpha.c.xi}: hyperparameter for learning \eqn{c^\xi}
#'      (argument is ignored when \code{learn.c.xi = FALSE}) (count component)
#'     \item \code{beta.c.tau}: hyperparameter for learning \eqn{c^\tau}
#'      (argument is ignored when \code{learn.c.tau = FALSE}) (count component)
#'     \item \code{beta.c.xi}: hyperparameter for learning \eqn{c^\xi}
#'      (argument is ignored when \code{learn.c.xi = FALSE}) (count component)
#'    \item \code{iota.c.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{c^\tau} (count component)
#'    \item \code{iota.c.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{c^\xi} (count component)
#'    \item \code{target.rate.c.tau}: desired acceptance rate when updating
#'     \eqn{c^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.c.tau = FALSE}) (count component)
#'    \item \code{target.rate.c.xi}: desired acceptance rate when updating
#'     \eqn{c^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.c.xi = FALSE}) (count component)
#'    \item \code{kappa.tau}: hyperparameter for learning \eqn{\tau^2_j}
#'     (double Gamma, original triple Gamma). Note that this is the global parameter.
#'     The component-specific parameters are learned automatically when triple
#'     Gamma prior is selected (count component)
#'    \item \code{kappa.xi}: hyperparameter for learning \eqn{\xi^2_j}
#'     (double Gamma, original triple Gamma). Note that this is the global parameter.
#'     The component-specific parameters are learned automatically when triple
#'     Gamma prior is selected (count component)
#'    \item \code{learn.kappa.tau}: if TRUE \eqn{\kappa^\tau} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.tau} used as starting value,
#'     if FALSE \eqn{\kappa^\tau = } \code{kappa.tau} (count component)
#'    \item \code{learn.kappa.xi}: if TRUE \eqn{\kappa^\xi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.xi} used as starting value,
#'     if FALSE \eqn{\kappa^\xi = } \code{kappa.xi} (count component)
#'    \item \code{d.tau}: hyperparameter for learning \eqn{\kappa^\tau}
#'     (only double Gamma) (count component)
#'    \item \code{e.tau}: hyperparameter for learning \eqn{\kappa^\tau}
#'     (double Gamma, original triple Gamma prior). Note that under triple Gamma shrinkage
#'     this parameter is only the starting value as it is learned during MCMC
#'     (count component)
#'    \item \code{d.xi}: hyperparameter for learning \eqn{\kappa^\xi}
#'     (only double Gamma) (count component)
#'    \item \code{e.xi}: hyperparameter for learning \eqn{\kappa^\xi}
#'     (double Gamma, original triple Gamma). Note that under triple Gamma shrinkage
#'     this parameter is only the starting value as it is learned during MCMC
#'      (count component)
#'    \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "rw-t0" (shrinkage prior starting at T = 0),
#'     "rw-t1" (shrinkage prior starting at T = 1) or "ind" (independence prior),
#'      (count component)
#'    \item \code{c}: prior parameter that scales the variance when using
#'     shrinkage prior (ignored when using independence prior) (count component)
#'    \item \code{B0}: prior variance on the regression parameters when using
#'     independence prior (ignored when using shrinkage prior) (count component)
#'    \item \code{TG}: Boolean indicating whether triple Gamma prior should be used;
#'     when \code{TG = FALSE} and \code{type = "rw-t0/1"} then the double Gamma
#'     prior is used; when \code{type = "ind"} the argument \code{TG} is ignored
#'      (count component)
#'    \item \code{TG.alternative}: Boolean indicating whether the alternative
#'     triple Gamma representation of Knaus and Frühwirth-Schnatter (2025) should
#'     be used. If \code{TG.alternative = FALSE}, the original triple Gamma
#'     representation based on Cadonna et al. (2020) will be used instead.
#'     Note that the hyperparameters \eqn{a} and \eqn{c} will not be sampled
#'     and have to be set to fixed values in the alternative representation.
#'     If \code{TG = FALSE}, this argument is ignored (count component)
#'   }
#' @param prior.load_zinb.count A list of arguments for estimating the parameters of the factor
#'  part of the count model. This argument is only used when \code{model = 'ZINB'} and
#'  otherwise ignored. The arguments are:
#'  \itemize{
#'   \item \code{a.phi}: hyperparameter for learning \eqn{\phi^2} (count component)
#'   \item \code{a.zeta}: hyperparameter for learning  \eqn{\zeta^2} (count component)
#'   \item \code{kappa.phi}: hyperparameter for learning \eqn{\phi^2} (count component)
#'   \item \code{kappa.zeta}: part of the rate parameter of the Gamma prior for \eqn{\zeta^2}
#'    (count component)
#'   \item \code{learn.kappa.phi}: if TRUE \eqn{\kappa^\phi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.phi} used as starting value,
#'     if FALSE \eqn{\kappa^\phi = } \code{kappa.phi} (count component)
#'   \item \code{learn.kappa.zeta}: if TRUE \eqn{\kappa^\zeta} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.zeta} used as starting value,
#'     if FALSE \eqn{\kappa^\zeta = } \code{kappa.zeta} (count component)
#'   \item \code{d.phi}: hyperparameter for learning \eqn{\kappa^\phi} (count component)
#'   \item \code{d.zeta}: hyperparameter for learning \eqn{\kappa^\zeta} (count component)
#'   \item \code{e.phi}: hyperparameter for learning \eqn{\kappa^\phi} (count component)
#'   \item \code{e.zeta}: hyperparameter for learning \eqn{\kappa^\zeta} (count component)
#'   \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "cps" (compound symmetric), "rw-t0" (shrinkage prior starting at T = 0),
#'     "rw-t1" (shrinkage prior starting at T = 1) or "ind" (independence prior)
#'      (count component)
#'   \item \code{c}: prior parameter that scales the variance when using shrinkage prior
#'    (ignored when using cps or independence prior) (count component)
#'   \item \code{L0}: prior variance on the factor loading when using either cps or independence
#'    prior (ignored when using shrinkage prior) (count component)
#'  }
#'  Note that for the factor model, the hyperparameters \code{a.phi} and
#'  \code{a.zeta} have to be fixed and are not sampled using Metropolis-Hastings.
#'  There is also no triple Gamma prior as only one factor is present
#' @param prior.reg_zinb.inflation A list of arguments for estimating the parameters of the regression
#'  part of the zero-inflation model. This argument is only used when \code{model = 'ZINB'} and
#'  otherwise ignored. The arguments are:
#'   \itemize{
#'    \item \code{a.tau}: hyperparameter for learning \eqn{\tau^2_j}
#'     (double Gamma, triple Gamma) (inflation component)
#'    \item \code{a.xi}: hyperparameter for learning \eqn{\xi^2_j}
#'     (double Gamma, triple Gamma) (inflation component)
#'    \item \code{learn.a.tau}: if TRUE \eqn{a^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.tau} used as starting value, if FALSE
#'     \eqn{a^\tau =} \code{a.tau} (not learned under alternative triple Gamma)
#'     (inflation component)
#'    \item \code{learn.a.xi}: if TRUE \eqn{a^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{a.xi} used as starting value, if FALSE
#'     \eqn{a^\xi =} \code{a.xi} (not learned under alternative triple Gamma)
#'      (inflation component)
#'     \item \code{alpha.a.tau}: hyperparameter for learning \eqn{a^\tau}
#'      (argument is ignored when \code{learn.a.tau = FALSE}) (inflation component)
#'     \item \code{alpha.a.xi}: hyperparameter for learning \eqn{a^\xi}
#'      (argument is ignored when \code{learn.a.xi = FALSE}) (inflation component)
#'     \item \code{beta.a.tau}: hyperparameter for learning \eqn{a^\tau}
#'      (argument is ignored when \code{learn.a.tau = FALSE}) (inflation component)
#'     \item \code{beta.a.xi}: hyperparameter for learning \eqn{a^\xi}
#'      (argument is ignored when \code{learn.a.xi = FALSE}) (inflation component)
#'    \item \code{iota.a.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\tau} (argument is ignored when
#'     \code{learn.a.tau = FALSE}) (inflation component)
#'    \item \code{iota.a.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{a^\xi} (argument is ignored when
#'     \code{learn.a.xi = FALSE}) (inflation component)
#'    \item \code{target.rate.a.tau}: desired acceptance rate when updating
#'     \eqn{a^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.tau = FALSE}) (inflation component)
#'    \item \code{target.rate.a.xi}: desired acceptance rate when updating
#'     \eqn{a^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.a.xi = FALSE}) (inflation component)
#'    \item \code{c.tau:} shape parameter in the Gamma prior for \eqn{\kappa^\tau_j}
#'     (only triple Gamma) (inflation component)
#'    \item \code{c.xi:} shape parameter in the Gamma prior for \eqn{\kappa^\xi_j}
#'     (only triple Gamma) (inflation component)
#'    \item \code{learn.c.tau}: if TRUE \eqn{c^\tau} is updated using Metropolis-Hastings
#'     with the value of argument \code{c.tau} used as starting value, if FALSE
#'     \eqn{c^\tau =} \code{c.tau} (not learned under alternative triple Gamma)
#'      (inflation component)
#'    \item \code{learn.c.xi}: if TRUE \eqn{c^\xi} is updated using Metropolis-Hastings
#'     with the value of argument \code{c.xi} used as starting value, if FALSE
#'     \eqn{c^\xi =} \code{c.xi} (not learned under alternative triple Gamma)
#'      (inflation component)
#'     \item \code{alpha.c.tau}: hyperparameter for learning \eqn{c^\tau}
#'      (argument is ignored when \code{learn.c.tau = FALSE}) (inflation component)
#'     \item \code{alpha.c.xi}: hyperparameter for learning \eqn{c^\xi}
#'      (argument is ignored when \code{learn.c.xi = FALSE}) (inflation component)
#'     \item \code{beta.c.tau}: hyperparameter for learning \eqn{c^\tau}
#'      (argument is ignored when \code{learn.c.tau = FALSE}) (inflation component)
#'     \item \code{beta.c.xi}: hyperparameter for learning \eqn{c^\xi}
#'      (argument is ignored when \code{learn.c.xi = FALSE}) (inflation component)
#'    \item \code{iota.c.tau}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{c^\tau} (inflation component)
#'    \item \code{iota.c.xi}: proposal standard deviation for Metropolis-Hastings
#'     updating of \eqn{c^\xi} (inflation component)
#'    \item \code{target.rate.c.tau}: desired acceptance rate when updating
#'     \eqn{c^\tau} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.c.tau = FALSE}) (inflation component)
#'    \item \code{target.rate.c.xi}: desired acceptance rate when updating
#'     \eqn{c^\xi} using Metropolis-Hastings (argument is ignored when
#'     \code{learn.c.xi = FALSE}) (inflation component)
#'    \item \code{kappa.tau}: hyperparameter for learning \eqn{\tau^2_j}
#'     (double Gamma, original triple Gamma). Note that this is the global parameter.
#'     The component-specific parameters are learned automatically when triple
#'     Gamma prior is selected (inflation component)
#'    \item \code{kappa.xi}: hyperparameter for learning \eqn{\xi^2_j}
#'     (double Gamma, original triple Gamma). Note that this is the global parameter.
#'     The component-specific parameters are learned automatically when triple
#'     Gamma prior is selected (inflation component)
#'    \item \code{learn.kappa.tau}: if TRUE \eqn{\kappa^\tau} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.tau} used as starting value,
#'     if FALSE \eqn{\kappa^\tau = } \code{kappa.tau} (inflation component)
#'    \item \code{learn.kappa.xi}: if TRUE \eqn{\kappa^\xi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.xi} used as starting value,
#'     if FALSE \eqn{\kappa^\xi = } \code{kappa.xi} (inflation component)
#'    \item \code{d.tau}: hyperparameter for learning \eqn{\kappa^\tau}
#'     (only double Gamma) (inflation component)
#'    \item \code{e.tau}: hyperparameter for learning \eqn{\kappa^\tau}
#'     (double Gamma, original triple Gamma prior). Note that under triple Gamma shrinkage
#'     this parameter is only the starting value as it is learned during MCMC
#'     (inflation component)
#'    \item \code{d.xi}: hyperparameter for learning \eqn{\kappa^\xi}
#'     (only double Gamma) (inflation component)
#'    \item \code{e.xi}: hyperparameter for learning \eqn{\kappa^\xi}
#'     (double Gamma, original triple Gamma). Note that under triple Gamma shrinkage
#'     this parameter is only the starting value as it is learned during MCMC
#'      (inflation component)
#'    \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "rw-t0" (shrinkage prior starting at T = 0),
#'     "rw-t1" (shrinkage prior starting at T = 1) or "ind" (independence prior),
#'      (inflation component)
#'    \item \code{c}: prior parameter that scales the variance when using
#'     shrinkage prior (ignored when using independence prior) (inflation component)
#'    \item \code{B0}: prior variance on the regression parameters when using
#'     independence prior (ignored when using shrinkage prior) (inflation component)
#'    \item \code{TG}: Boolean indicating whether triple Gamma prior should be used;
#'     when \code{TG = FALSE} and \code{type = "rw-t0/1"} then the double Gamma
#'     prior is used; when \code{type = "ind"} the argument \code{TG} is ignored
#'      (inflation component)
#'    \item \code{TG.alternative}: Boolean indicating whether the alternative
#'     triple Gamma representation of Knaus and Frühwirth-Schnatter (2025) should
#'     be used. If \code{TG.alternative = FALSE}, the original triple Gamma
#'     representation based on Cadonna et al. (2020) will be used instead.
#'     Note that the hyperparameters \eqn{a} and \eqn{c} will not be sampled
#'     and have to be set to fixed values in the alternative representation.
#'     If \code{TG = FALSE}, this argument is ignored (inflation component)
#'   }
#' @param prior.load_zinb.inflation A list of arguments for estimating the parameters of the factor
#'  part of the zero-inflation model. This argument is only used when \code{model = 'ZINB'} and
#'  otherwise ignored. The arguments are:
#'  \itemize{
#'   \item \code{a.phi}: hyperparameter for learning \eqn{\phi^2} (inflation component)
#'   \item \code{a.zeta}: hyperparameter for learning  \eqn{\zeta^2} (inflation component)
#'   \item \code{kappa.phi}: hyperparameter for learning \eqn{\phi^2} (inflation component)
#'   \item \code{kappa.zeta}: part of the rate parameter of the Gamma prior for \eqn{\zeta^2}
#'    (inflation component)
#'   \item \code{learn.kappa.phi}: if TRUE \eqn{\kappa^\phi} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.phi} used as starting value,
#'     if FALSE \eqn{\kappa^\phi = } \code{kappa.phi} (inflation component)
#'   \item \code{learn.kappa.zeta}: if TRUE \eqn{\kappa^\zeta} is sampled in a
#'     Gibbs-step with the value of argument \code{kappa.zeta} used as starting value,
#'     if FALSE \eqn{\kappa^\zeta = } \code{kappa.zeta} (inflation component)
#'   \item \code{d.phi}: hyperparameter for learning \eqn{\kappa^\phi} (inflation component)
#'   \item \code{d.zeta}: hyperparameter for learning \eqn{\kappa^\zeta} (inflation component)
#'   \item \code{e.phi}: hyperparameter for learning \eqn{\kappa^\phi} (inflation component)
#'   \item \code{e.zeta}: hyperparameter for learning \eqn{\kappa^\zeta} (inflation component)
#'   \item \code{type}: the type of prior you want on your regression effects;
#'     this argument is either "cps" (compound symmetric), "rw-t0" (shrinkage prior starting at T = 0),
#'     "rw-t1" (shrinkage prior starting at T = 1) or "ind" (independence prior)
#'      (inflation component)
#'   \item \code{c}: prior parameter that scales the variance when using shrinkage prior
#'    (ignored when using cps or independence prior) (inflation component)
#'   \item \code{L0}: prior variance on the factor loading when using either cps or independence
#'    prior (ignored when using shrinkage prior) (inflation component)
#'  }
#'  Note that for the factor model, the hyperparameters \code{a.phi} and
#'  \code{a.zeta} have to be fixed and are not sampled using Metropolis-Hastings.
#'  There is also no triple Gamma prior as only one factor is present
#' @param mcmc.opt a list containing information on the overall sampler.
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
#'  parameter \eqn{r} in the Negative Binomial and Zero-Inflated Negative Binomial
#'  model using univariate Slice sampling.
#'  For other response distributions, this is ignored. The arguments are:
#'  \itemize{
#'    \item \code{alpha.r}: shape parameter of Gamma proposal
#'    \item \code{beta.r}: rate parameter of Gamma proposal
#'    \item \code{expansion.steps}: number of steps in stepping-out phase
#'    \item \code{width}: width of the slice interval (on the log-scale of \eqn{r})
#'    \item \code{p.overrelax}: probability of performing an overrelaxation step;
#'     performing overrelaxation might increase sampling efficiency; when overrelaxation
#'     should not be used, set \code{p.overrelax = 0}
#'    \item \code{accuracy.overrelax}: accuracy in overrelaxation phase
#'  }
#'  For more information on overrelaxation and Slice sampling in general,
#'  see the original paper on Slice sampling by Neal (2003)
#' @param HPD.coverage coverage probability of highest posterior density intervals
#'  (default yields 95 percent coverage)
#' @param R.WAIC number of replications for computing the marginal WAIC, where
#'  the latent factors are integrated out
#' @param posterior.predictive.matrix if TRUE (= default) the posterior predictive distribution
#'  based on the training data is computed and returned as a list object by the function.
#'  Setting it to FALSE is usually only done for saving memory
#' @param random.effects if TRUE (= default) a factor model is included for estimating
#'  random effects, if FALSE the model does not contain random effects and, consequently,
#'  priors on the parameters of the factor model are ignored
#' @param progress.bar if TRUE a progress bar is displayed, if FALSE the progress bar is omitted
#'
#' @description
#' This is the main function for fitting a flexible Bayesian panel data model
#'  in the time-varying parameter framework. By using shrinkage priors, it is
#'  possible to identify whether an effect is time-varying, time-invariant or zero.
#'  This function works for Gaussian, binary and (zero-inflated) Negative Binomial response data.
#'
#' @details
#'  This function fits a Bayesian time-varying parameter panel data model to a
#'  longitudinal response \eqn{\text{y}_{it}} for \eqn{i \in \{1,\dots,n\}} subjects that are observed
#'  at \eqn{t \in \{1,\dots,T\}}
#'  time points. The model is expressed in its non-centered form
#'  (see Frühwirth-Schnatter and Wagner (2010) for details on this parameterization)
#'  as this makes it easier to place shrinkage priors on the model parameters.
#'  By using the non-centered parameterization, the response is linked to the
#'  following linear predictor
#'  \deqn{\eta_{it} = \mathbf{x}_{it}^\top \boldsymbol{\beta} +
#'  \mathbf{x}_{it}^\top \boldsymbol{\Theta} \boldsymbol{\tilde{\beta}}_t +
#'   f_i\lambda + f_i \psi \tilde{\lambda}_t,}
#'  where \eqn{\mathbf{x}_{it}} is a \eqn{d}-dimensional column vector of covariate
#'  values for subject \eqn{i} at time \eqn{t}, \eqn{\boldsymbol{\beta}} is a \eqn{d}-dimensional
#'  column vector of fixed effects (including the intercept as first parameter),
#'  \eqn{\boldsymbol{\Theta} = \text{diag}(\theta_1,\dots,\theta_d)} is a
#'  diagonal matrix, where larger main diagonal elements indicate stronger variation
#'  of the regression effects over time, \eqn{\boldsymbol{\tilde{\beta}}_t}
#'  is a \eqn{d}-dimensional state vector that follows a Gaussian random walk,
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
#'  prevents the model from overfitting, i.e., it is reasonable to assume that not
#'  every covariate has a time-varying effect and without proper regularization
#'  estimates are likely less stable.
#'  The priors on the parameters are specified as Normal-Gamma (or double Gamma)
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
#'
#'  We follow Bitto and Frühwirth-Schnatter (2019) and assign the following hyperpriors
#'  on the parameters of the double Gamma prior
#'  \deqn{
#'  \begin{aligned}
#'   \kappa^\xi|d^\xi,e^\xi& \sim \mathcal{G}(d^\xi,e^\xi), \quad  a^\xi|\alpha_{a^\xi},\beta_{a^\xi} \sim  \mathcal{G}(\alpha_{a^\xi}, \alpha_{a^\xi} \beta_{a^\xi}), \\
#'   \kappa^\tau|d^\tau,e^\tau& \sim  \mathcal{G}(d^\tau,e^\tau), \quad a^\tau|\alpha_{a^\tau},\beta_{a^\tau} \sim  \mathcal{G}(\alpha_{a^\tau}, \alpha_{a^\tau} \beta_{a^\tau}), \\
#'   \kappa^\zeta|d^\zeta,e^\zeta& \sim \mathcal{G}(d^\zeta,e^\zeta), \\
#'   \kappa^\phi|d^\phi,e^\phi& \sim \mathcal{G}(d^\phi,e^\phi).
#'  \end{aligned}
#' }
#'  Note that for \eqn{a \le 1} the priors are valid shrinkage priors. The
#'  Bayesian Lasso is a special case of the Normal-Gamma prior with \eqn{a = 1}.
#'  In our simulations we have achieved good results by using the Bayesian Lasso
#'  for the hyperparameters in the factor model (\eqn{a^\phi, a^\zeta}) and by sampling
#'  the hyperparameters in the regression model (\eqn{a^\tau, a^\xi}) using
#'  Metropolis-Hastings updates. For inference, we have adapted the MCMC sampler
#'  presented in Bitto and Frühwirth-Schnatter (2019) and implemented in the \code{R} package
#'  \code{shrinkTVP} (Knaus et al., 2021).
#'
#'  We also consider an alternative and more general shrinkage prior - the triple Gamma
#'  prior - which encompasses the double Gamma prior as a special case (Cadonna et al., 2020).
#'  The triple Gamma is a Normal-Gamma-Gamma prior on the scale and regression parameters, i.e.,
#'  \deqn{
#'  \begin{aligned}
#'    \theta_j|\xi_j^2 &\sim \mathcal{N}(0,\xi^2_j),
#'    &\quad \xi_j^2|a^\xi, \kappa^\xi_j &\sim \mathcal{G}\left(a^\xi, \frac{a^\xi \kappa_j^\xi}{2}\right),
#'   &\quad \kappa_j^\xi|c^\xi,\kappa^\xi
#'     \sim \mathcal{G}\left(c^\xi, \frac{c^\xi}{\kappa^\xi}\right) &\quad j = \{1, \dots, d\}, \\
#'    \beta_j|\tau_j^2 &\sim \mathcal{N}(0,\tau_j^2),
#'    &\quad \tau_j^2|a^\tau, \kappa^\tau_j &\sim \mathcal{G}\left(a^\tau, \frac{a^\tau \kappa_j^\tau}{2}\right),
#'    &\quad \kappa_j^\tau|c^\tau,\kappa^\tau
#'     \sim \mathcal{G}\left(c^\tau, \frac{c^\tau}{\kappa^\tau}\right)
#'    &\quad j = \{1, \dots, d\}. \\
#'  \end{aligned}
#'  }
#'  The additional prior-layer allows for covariate specific hyperparameters.
#'  As the factor part of the model contains only one factor, the triple Gamma
#'  prior only makes sense for the regression part of the model.
#'
#'  For building the MCMC sampler, we follow Cadonna et al. (2020) and consider
#'  the following, alternative representation of the triple Gamma prior as the
#'  basis for inference
#'  \deqn{
#'    \begin{aligned}
#'    \theta_j|\overset{\vee}{\xi^2_j},\overset{\vee}{\kappa_j^\xi},\phi^\xi &\sim
#'     \mathcal{N}(0,\phi^\xi\overset{\vee}{\xi_j^2}/\overset{\vee}{\kappa_j^\xi}),
#'     \quad \overset{\vee}{\xi_j^2}|a^\xi \sim \mathcal{G}(a^\xi,1), \quad
#'     \overset{\vee}{\kappa_j^\xi}|c^\xi \sim \mathcal{G}(c^\xi,1), \\
#'     \beta_j|\overset{\vee}{\tau^2_j},\overset{\vee}{\kappa_j^\tau},\phi^\tau &\sim
#'     \mathcal{N}(0,\phi^\tau\overset{\vee}{\tau_j^2}/\overset{\vee}{\kappa_j^\tau}),
#'     \quad \overset{\vee}{\tau_j^2}|a^\tau \sim \mathcal{G}(a^\tau,1), \quad
#'     \overset{\vee}{\kappa_j^\tau}|c^\tau \sim \mathcal{G}(c^\tau,1),
#'    \end{aligned}
#'  }
#'  where setting \eqn{\tau^2_j = \phi^\tau\overset{\vee}{\tau_j^2}/\overset{\vee}{\kappa_j^\tau}}
#'  and \eqn{\xi^2_j = \phi^\xi\overset{\vee}{\xi_j^2}/\overset{\vee}{\kappa_j^\xi}}
#'  makes the connection to the original triple Gamma representation more obvious.
#'  Here, \eqn{\phi^\tau = 2c^\tau/(a^\tau \kappa^\tau)} and
#'  \eqn{\phi^\xi = 2c^\xi/(a^\xi \kappa^\xi)}.
#'
#'  Under triple Gamma shrinkage, the priors for \eqn{a^\xi,a^\tau,c^\xi,c^\tau}
#'  are selected such that their support is restricted to \eqn{(0,0.5)} as the
#'  triple Gamma prior then has a pole at the origin.
#'   Cadonna et al. (2020) then consider the following Beta-priors
#'  \deqn{
#'    \begin{aligned}
#'     2a^\xi|\alpha_{a^\xi},\beta_{a^\xi} &\sim \mathcal{B}(\alpha_{a^\xi},\beta_{a^\xi}), \quad 2c^\xi|\alpha_{c^\xi}, \beta_{c^\xi} \sim \mathcal{B}(\alpha_{c^\xi},\beta_{c^\xi}) \\
#'     2a^\tau|\alpha_{a^\tau}, \beta_{a^\tau} &\sim \mathcal{B}(\alpha_{a^\tau},\beta_{a^\tau}), \quad 2c^\tau|\alpha_{c^\tau}, \beta_{c^\tau} \sim \mathcal{B}(\alpha_{c^\tau},\beta_{c^\tau})
#'    \end{aligned}
#'  }
#'  Moreover, the global shrinkage parameters are equipped with the following F-priors
#'  \deqn{
#'    \frac{\kappa^\xi}{2}|a^\xi,c^\xi \sim \mathcal{F}(2a^\xi,2c^\xi), \quad
#'    \frac{\kappa^\tau}{2}|a^\tau,c^\tau \sim \mathcal{F}(2a^\tau,2c^\tau).
#'  }
#'  Important: To be consistent with the double Gamma prior, the parameters
#'  \eqn{\xi^2_j} and \eqn{\tau^2_j} are returned by the function, whereas
#'  \eqn{\overset{\vee}{\xi_j^2}} and \eqn{\overset{\vee}{\tau_j^2}} are omitted
#'  from the output and only sampled internally during MCMC. However, they can
#'  be easily computed based on the above definitions.
#'
#'  It should be stated that not only the double Gamma prior and the Bayesian
#'  Lasso prior are special cases of the triple Gamma prior, but also other well-known
#'  shrinkage priors such as the Strawderman-Berger prior and the Horseshoe prior.
#'  For obtaining the Strawderman-Berger prior, set \eqn{a=0.5, c = 1, \kappa = 4}.
#'  For obtaining the Horseshoe prior, set \eqn{a =0.5, c = 0.5}.
#'  The \eqn{\kappa} parameters can then either be learned (hierachical Horseshoe)
#'  or held fixed (non-hierarchical Horseshoe) (see Cadonna et al., 2020).
#'
#'  Finally, there exists an alternative representation of the triple Gamma prior,
#'  which was proposed by Knaus and Frühwirth-Schnatter (2025) and is given by
#'  \deqn{
#'   \begin{aligned}
#'    \theta_j|\xi^2_j &\sim \mathcal{N}(0,\xi^2_j), \quad \xi^2_j|\chi^\xi_j,c^\xi \sim
#'    \mathcal{G}^{-1}(c^\xi,\chi^\xi_j), \quad \chi_j^\xi|a^\xi,c^\xi \sim
#'    \mathcal{G}\Biggl(a^\xi, \frac{a^\xi}{c^\xi}\Biggr), \\
#'    \beta_j|\tau^2_j &\sim \mathcal{N}(0,\tau^2_j), \quad \tau^2_j|\chi^\tau_j,c^\tau \sim
#'    \mathcal{G}^{-1}(c^\tau,\chi^\tau_j), \quad \chi_j^\tau|a^\tau,c^\tau \sim
#'    \mathcal{G}\Biggl(a^\tau, \frac{a^\tau}{c^\tau}\Biggr).
#'   \end{aligned}
#'  }
#'  This new representation of the Triple Gamma shrinkage prior fascilitates posterior inference
#'  as the parameters are sampled from standard distributions.
#'  Following Knaus and Frühwirth-Schnatter (2025), the hyperparameters \eqn{a} and \eqn{c}
#'  are assumed to be fixed by the researcher.
#'
#'  The function \code{panelTVP} can handle the following popular regression models
#'   \itemize{
#'     \item Gaussian model (Normal outcomes)
#'     \item Logit model (binary outcomes)
#'     \item Probit model (binary outcomes)
#'     \item Negative Binomial model (overdispersed count outcomes)
#'     \item Zero-Inflated Negative Binomial model (both zero-inflated and overdispersed count outcomes)
#'   }
#'
#'  With a Gaussian (Normal) response, the model is given as follows
#'  \deqn{\text{y}_{it} = \eta_{it} + \varepsilon_{it}, \quad \varepsilon_{it} \sim \mathcal{N}(0,\sigma^2),}
#'  where we place a hierarchical prior on the error variance
#'  \deqn{\sigma^2|c_0,C_0 \sim \mathcal{G}^{-1}(c_0,C_0), \quad C_0|g_0,G_0 \sim \mathcal{G}(g_0,G_0).}
#'
#'  In the binary Logit model, the conditional probability that \eqn{\text{y}_{it} = 1} is modelled
#'  via the Logit-link as
#'  \deqn{\mathbb{P}(\text{y}_{it} = 1|\eta_{it}) = \frac{\exp(\eta_{it})}{1+\exp(\eta_{it})},   }
#'  whereas in the binary Probit model, this probability is modelled
#'  via the Probit-link as
#'  \deqn{\mathbb{P}(\text{y}_{it} = 1|\eta_{it}) = \Phi(\eta_{it}),}
#'  where \eqn{\Phi(\cdot)} denotes the standard Normal cumulative distribution function.
#'
#'  For modelling overdispersed count data, we assume that the response is a realization
#'  of a Negative Binomial distribution with pmf given by
#'  \deqn{p(\text{y}_{it}|r,\eta_{it}) = \frac{\Gamma(\text{y}_{it}+r)}{\Gamma(r)\text{y}_{it}!}(1-q_{it})^r q_{it}^{\text{y}_{it}},
#'   \quad q_{it} = \frac{\exp(\eta_{it})}{1+\exp(\eta_{it})},}
#'  where we place a prior on the constant dispersion parameter
#'   \deqn{r|\alpha^r,\beta^r \sim \mathcal{G}(\alpha^r,\beta^r).}
#'
#'  Inference in the Logit and Negative Binomial model is based on data augmentation using
#'  Pólya-Gamma random variables (see Pillow and Scott, 2012; Polson et al., 2013).
#'
#'  For modelling both zero-inflated and overdispersed count
#'  data, we assume that there are two latent classes of zeros to account for the excess zeros:
#'  structural and at-risk zeros. A structural zero belongs to an observation that is not at risk of
#'  experiencing the event, whereas an at-risk zero belongs to an observation that is
#'  at-risk of experiencing the event but has for some reason a zero recorded.
#'  It is therefore assumed that the outcome in the Zero-Inflated Negative Binomial model
#'   \eqn{\text{y}_{it}} is a realization of a
#'  mixture distribution with a point mass at zero (for the structural zeros) and
#'  a standard Negative Binomial count model for observations that are at-risk,
#'  i.e., this includes both at-risk zeros and positive counts. The ZINB model
#'  can, thus, be stated as (Neelon, 2019)
#'  \deqn{ \text{y}_{it}|r,\mu_{it},w_{it} \sim (1-\pi_{it}) \cdot \mathbb{I}_{(w_{it} = 0)} + \pi_{it} \cdot \mathcal{NB}(\mu_{it},r) \mathbb{I}_{(w_{it}=1)},}
#'  where \eqn{w_{it}} is a latent at-risk indicator such that with probability
#'  \eqn{1-\pi_{it}, w_{it} = 0} and with probability \eqn{\pi_{it}, w_{it} = 1},
#'  and \eqn{\mu_{it}} is the mean of the Negative Binomial distribution.
#'
#'  There are two separate linear predictors in the ZINB model. The first
#'  linear predictor \eqn{\eta_{it}^{\text{logit}}} is related to the zero-inflation
#'  part of the model through
#'  \deqn{\pi_{it} = \frac{\exp(\eta_{it}^{\text{logit}})}{1+\exp(\eta_{it}^{\text{logit}})},}
#'  whereas the second linear predictor \eqn{\eta_{it}^{\text{nb}}} is related to
#'  the count part of the model through
#'  \deqn{q_{it} = \frac{\exp(\eta_{it}^{\text{nb}})}{1+\exp(\eta_{it}^{\text{nb}})}.}
#'  Moreover, the mean of the count component is given as
#'  \deqn{\mu_{it} = r \exp(\eta_{it}^{\text{nb}}).}
#'  In every iteration of the MCMC sampler, the latent at-risk indicators \eqn{w_{it}} are updated
#'  and the parameters in the count model are estimated by using only the observations
#'  that are currently in the risk set. Note that for \eqn{\text{y}_{it} > 0} the at-risk indicators
#'  are fixed at \eqn{w_{it} = 1} in every iteration.
#'
#'  If the response variable contains missing values, data augmentation is used to
#'  impute those values. This works for all models regardless of the distributional assumption.
#'
#' @returns The function returns an object of class \code{panelTVP.Gaussian},
#'  \code{panelTVP.Probit}, \code{panelTVP.Logit}, \code{panelTVP.NegBin} or
#'  \code{panelTVP.ZINB}
#'  depending on whether a Gaussian, Probit, Logit, Negative Binomial or
#'  Zero-Inflated Negative Binomial model was
#'  fitted. When either a Gaussian, Probit, Logit or Negative Binomial model was fitted,
#'  the returned object contains a list of the following elements:
#'  \describe{
#'    \item{learning.settings}{information on which parameters have been learned}
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
#'    \item{fmean}{posterior means of factor scores (only when random effects are estimated)}
#'    \item{model}{the fitted model}
#'    \item{acceptance.rates}{the achieved acceptance rates when using Metropolis-Hastings}
#'    \item{HPD.coverage}{coverage probability of HPD intervals (based on input)}
#'    \item{runtime}{total runtime of the sampler (measured in seconds)}
#'    \item{WAIC}{the Widely Applicable Information Criterion for model comparison
#'     (note that the WAIC is only computed from the actually observed data,
#'     i.e., missing data are fully ignored when computing WAIC)}
#'    \item{posterior.predictive}{the \eqn{Tn \times M} matrix that contains the posterior
#'     predictive distribution for each observation (rows) and each MCMC draw (columns)}
#'    \item{mcmc.settings}{details on the MCMC sampler}
#'    \item{variable.codes}{information on which covariate is associated with which parameter}
#'  }
#'  When modelling a Zero-Inflated Negative Binomial response,
#'  the returned object contains a list of the following elements:
#'  \describe{
#'    \item{learning.settings_logit}{information on which parameters have been learned
#'       for the zero-inflation component of the model}
#'    \item{learning.settings_nb}{information on which parameters have been learned
#'       for the count component of the model}
#'    \item{data}{the data used for fitting the model and additional context information
#'    derived from the data}
#'    \item{Y}{the \eqn{Tn \times M} response data matrix of every iteration of the chain,
#'    i.e., this matrix is only included in the output when missing response data
#'    were present and imputed via data augmentation. Each row of \code{Y} contains
#'    the sampled values of a specific observation. For observed data, the corresponding
#'    rows are essentially replicates of the same value. For missing data, the
#'    corresponding rows contain the imputed values of every iteration.}
#'    \item{mcmc_logit}{Markov Chains for every parameter except for the factor scores
#'     (to save memory) for the zero-inflation component of the model}
#'    \item{mcmc_nb}{Markov Chains for every parameter except for the factor scores
#'     (to save memory) for the count component of the model}
#'    \item{posterior_logit}{preliminary summary of posterior results for the
#'     zero-inflation component of the model}
#'    \item{posterior_nb}{preliminary summary of posterior results for the
#'     count component of the model}
#'    \item{fmean_logit}{posterior means of factor scores for the
#'      zero-inflation component of the model (only when random effects are estimated)}
#'    \item{fmean_nb}{posterior means of factor scores for the
#'      count component of the model (only when random effects are estimated)}
#'    \item{model}{the fitted model}
#'    \item{acceptance.rates}{the achieved acceptance rates when using
#'       Metropolis-Hastings for both components of the model}
#'    \item{HPD.coverage}{coverage probability of HPD intervals (based on input)}
#'    \item{runtime}{total runtime of the sampler (measured in seconds)}
#'    \item{WAIC}{the Widely Applicable Information Criterion for model comparison
#'     (note that the WAIC is only computed from the actually observed data,
#'     i.e., missing data are fully ignored when computing WAIC)}
#'    \item{posterior.predictive}{the \eqn{Tn \times M} matrix that contains the posterior
#'     predictive distribution for each observation (rows) and each MCMC draw (columns)}
#'    \item{mcmc.settings}{details on the MCMC sampler}
#'    \item{variable.codes_logit}{information on which covariate is associated with which parameter
#'       for the zero-inflation component of the model}
#'    \item{variable.codes_nb}{information on which covariate is associated with which parameter
#'       for the count component of the model}
#'  }
#' @export
#'
#' @author Roman Pfeiler, Helga Wagner
#'
#' @examples
#' # Example 1: Gaussian Response, learn all possible parameters (double Gamma prior)
#' sim.gaussian <- sim_panelTVP(n = 100,
#'                              Tmax = 4,
#'                              beta = c(4,1,0,0),
#'                              theta = c(1,0.5,0,0),
#'                              lambda = 1,
#'                              psi = 0.2,
#'                              model = "Gaussian",
#'                              sigma2 = 0.7)
#' res.gaussian1 <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#' # Example 2: Gaussian Response, hierarchical Bayesian Lasso (double Gamma prior)
#' prior.reg <- utils::modifyList(as.list(rlang::fn_fmls(panelTVP)$prior.reg),
#'                                list(a.tau = 1, a.xi = 1,
#'                                     learn.a.tau = FALSE, learn.a.xi = FALSE))
#' prior.load <- utils::modifyList(as.list(rlang::fn_fmls(panelTVP)$prior.load),
#'                                list(a.phi = 1, a.zeta = 1,
#'                                     learn.a.phi = FALSE, learn.a.zeta = FALSE))
#' res.gaussian2 <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          prior.reg = prior.reg,
#'                          prior.load = prior.load,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#'
#' # Example 3: Gaussian Response, learn all possible parameters (original triple Gamma prior)
#' prior.reg <- utils::modifyList(as.list(rlang::fn_fmls(panelTVP)$prior.reg), list(TG = TRUE))
#' res.gaussian3 <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          prior.reg = prior.reg,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#'
#' # Example 4: Gaussian Response, hierarchical Horseshoe prior (original triple Gamma prior)
#' prior.reg <- utils::modifyList(as.list(rlang::fn_fmls(panelTVP)$prior.reg),
#'                                list(TG = TRUE, a.tau = 0.5, a.xi = 0.5,
#'                                     learn.a.tau = FALSE, learn.a.xi = FALSE,
#'                                     c.tau = 0.5, c.xi = 0.5,
#'                                     learn.c.tau = FALSE, learn.c.xi = FALSE))
#' res.gaussian4 <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          prior.reg = prior.reg,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#'
#' # Example 5: Gaussian Response, non-hierarchical Horseshoe prior (original triple Gamma prior)
#' prior.reg <- utils::modifyList(as.list(rlang::fn_fmls(panelTVP)$prior.reg),
#'                                list(TG = TRUE, a.tau = 0.5, a.xi = 0.5,
#'                                     learn.a.tau = FALSE, learn.a.xi = FALSE,
#'                                     c.tau = 0.5, c.xi = 0.5,
#'                                     learn.c.tau = FALSE, learn.c.xi = FALSE,
#'                                     learn.kappa.tau = FALSE, learn.kappa.xi = FALSE))
#' res.gaussian5 <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          prior.reg = prior.reg,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#'
#' # Example 6: Gaussian Response, Strawderman-Berger prior (original triple Gamma prior)
#' prior.reg <- utils::modifyList(as.list(rlang::fn_fmls(panelTVP)$prior.reg),
#'                                list(TG = TRUE, a.tau = 0.5, a.xi = 0.5,
#'                                     learn.a.tau = FALSE, learn.a.xi = FALSE,
#'                                     c.tau = 1, c.xi = 1,
#'                                     learn.c.tau = FALSE, learn.c.xi = FALSE,
#'                                     kappa.tau = 4, kappa.xi = 4,
#'                                     learn.kappa.tau = FALSE, learn.kappa.xi = FALSE))
#' res.gaussian6 <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          prior.reg = prior.reg,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#'
#' # Example 7: Gaussian Response, alternative triple Gamma prior
#' prior.reg <- utils::modifyList(as.list(rlang::fn_fmls(panelTVP)$prior.reg),
#'                                list(TG = TRUE, TG.alternative = TRUE))
#' res.gaussian7 <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          prior.reg = prior.reg,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#'
#' # Example 8: Logit Response, learn all possible parameters (double Gamma prior)
#' sim.logit <- sim_panelTVP(n = 100,
#'                           Tmax = 4,
#'                           beta = c(1,0.5,0,0),
#'                           theta = c(0.8,0.5,0,0),
#'                           lambda = 1,
#'                           psi = 0.2,
#'                           model = "Logit")
#' res.logit <- panelTVP(y ~ W1 + W2 + W3,
#'                       data = sim.logit$observed,
#'                       id = sim.logit$observed$id,
#'                       t = sim.logit$observed$t,
#'                       mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                       model = "Logit")
#'
#' # Example 9: ZINB, learn all possible parameters (double Gamma prior)
#' sim.zinb <- sim_panelTVP(n = 100,
#'                          Tmax = 4,
#'                          beta_zinb.count = c(0.5,-0.7,0,0),
#'                          theta_zinb.count = c(0.05,0.5,0,0),
#'                          lambda_zinb.count = 0.5,
#'                          psi_zinb.count = 0.02,
#'                          beta_zinb.inflation = c(-1,0.6,0,0),
#'                          theta_zinb.inflation = c(0,1,0,0),
#'                          lambda_zinb.inflation = 0.7,
#'                          psi_zinb.inflation = 0,
#'                          r = 2,
#'                          model = "ZINB")
#' res.zinb <- panelTVP(y ~ W1.nb + W2.nb + W3.nb | W1.logit + W2.logit + W3.logit,
#'                      data = sim.zinb$observed,
#'                      id = sim.zinb$observed$id,
#'                      t = sim.zinb$observed$t,
#'                      mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                      model = "ZINB")
#'
#' @references
#'
#'  Bitto, A. and Frühwirth-Schnatter, S. (2019). Achieving Shrinkage in a
#'  Time-Varying Parameter Model Framework. In: Journal of Econometrics, 210,
#'  75-97.
#'
#'  Cadonna, A., Frühwirth-Schnatter, S. and Knaus, P. (2020). Triple the Gamma -
#'  A Unifying Shrinkage Prior for Variance and Variable Selection in Sparse
#'  State Space and TVP Models. In: Econometrics, 8, 1-36.
#'
#'  Frühwirth-Schnatter, S. and Wagner, H. (2010). Stochastic Model Specification
#'  Search for Gaussian and Partially Non-Gaussian State Space Models. Journal
#'  of Econometrics, 154, 85-100.
#'
#'  Knaus, P., Bitto-Nemling, A., Cadonna, A. and Frühwirth-Schnatter, S. (2021).
#'  Shrinkage in the Time-Varying Parameter Model Framework using the R Package
#'  \code{shrinkTVP}. In: Journal of Statistical Software, 100, 1-32.
#'
#'  Knaus, P. Frühwirth-Schnatter, S. (2025). The Dynamic Triple Gamma as a Shrinkage
#'  Process for Time-Varying Parameter Models, arXiv preprint arXiv:2312.10487v2.
#'
#'  Neal, R.M. (2003). Slice sampling. In: The Annals of Statistics, 31, 705-767.
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
#' @import stats
panelTVP <- function(formula = NULL,
                     data = NULL,
                     id = NULL,
                     t = NULL,
                     model = NULL,
                     prior.reg = list(
                       a.tau = 0.1, a.xi = 0.1,
                       learn.a.tau = TRUE, learn.a.xi = TRUE,
                       alpha.a.tau = 2, alpha.a.xi = 2,
                       beta.a.tau = 1, beta.a.xi = 1,
                       iota.a.tau = 1, iota.a.xi = 1,
                       target.rate.a.tau = 0.44, target.rate.a.xi = 0.44,
                       c.tau = 0.1, c.xi = 0.1,
                       learn.c.tau = TRUE, learn.c.xi = TRUE,
                       alpha.c.tau = 2, alpha.c.xi = 2,
                       beta.c.tau = 1, beta.c.xi = 1,
                       iota.c.tau = 1, iota.c.xi = 1,
                       target.rate.c.tau = 0.44, target.rate.c.xi = 0.44,
                       kappa.tau = 10, kappa.xi = 10,
                       learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                       d.tau = 0.001, d.xi = 0.001,
                       e.tau = 0.001, e.xi = 0.001,
                       type = "rw-t1", c = 1, B0 = 1, TG = FALSE, TG.alternative = FALSE
                     ),
                     prior.var = list(
                       learn.C0.hyp = list(g0 = 5, G0 = 3.333333), c0 = 2.5
                     ),
                     prior.load = list(
                       a.phi = 0.1, a.zeta = 0.1,
                       kappa.phi = 10, kappa.zeta = 10,
                       learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                       d.phi = 0.001, d.zeta = 0.001,
                       e.phi = 0.001, e.zeta = 0.001,
                       type = "rw-t1", c = 1, L0 = 1
                     ),
                     prior.reg_zinb.count = list(
                       a.tau = 0.1, a.xi = 0.1,
                       learn.a.tau = TRUE, learn.a.xi = TRUE,
                       alpha.a.tau = 2, alpha.a.xi = 2,
                       beta.a.tau = 1, beta.a.xi = 1,
                       iota.a.tau = 1, iota.a.xi = 1,
                       target.rate.a.tau = 0.44, target.rate.a.xi = 0.44,
                       c.tau = 0.1, c.xi = 0.1,
                       learn.c.tau = TRUE, learn.c.xi = TRUE,
                       alpha.c.tau = 2, alpha.c.xi = 2,
                       beta.c.tau = 1, beta.c.xi = 1,
                       iota.c.tau = 1, iota.c.xi = 1,
                       target.rate.c.tau = 0.44, target.rate.c.xi = 0.44,
                       kappa.tau = 10, kappa.xi = 10,
                       learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                       d.tau = 0.001, d.xi = 0.001,
                       e.tau = 0.001, e.xi = 0.001,
                       type = "rw-t1", c = 1, B0 = 1, TG = FALSE, TG.alternative = FALSE
                     ),
                     prior.load_zinb.count = list(
                       a.phi = 0.1, a.zeta = 0.1,
                       kappa.phi = 10, kappa.zeta = 10,
                       learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                       d.phi = 0.001, d.zeta = 0.001,
                       e.phi = 0.001, e.zeta = 0.001,
                       type = "rw-t1", c = 1, L0 = 1
                     ),
                     prior.reg_zinb.inflation = list(
                       a.tau = 0.1, a.xi = 0.1,
                       learn.a.tau = TRUE, learn.a.xi = TRUE,
                       alpha.a.tau = 2, alpha.a.xi = 2,
                       beta.a.tau = 1, beta.a.xi = 1,
                       iota.a.tau = 1, iota.a.xi = 1,
                       target.rate.a.tau = 0.44, target.rate.a.xi = 0.44,
                       c.tau = 0.1, c.xi = 0.1,
                       learn.c.tau = TRUE, learn.c.xi = TRUE,
                       alpha.c.tau = 2, alpha.c.xi = 2,
                       beta.c.tau = 1, beta.c.xi = 1,
                       iota.c.tau = 1, iota.c.xi = 1,
                       target.rate.c.tau = 0.44, target.rate.c.xi = 0.44,
                       kappa.tau = 10, kappa.xi = 10,
                       learn.kappa.tau = TRUE, learn.kappa.xi = TRUE,
                       d.tau = 0.001, d.xi = 0.001,
                       e.tau = 0.001, e.xi = 0.001,
                       type = "rw-t1", c = 1, B0 = 1, TG = FALSE, TG.alternative = FALSE
                     ),
                     prior.load_zinb.inflation = list(
                       a.phi = 0.1, a.zeta = 0.1,
                       kappa.phi = 10, kappa.zeta = 10,
                       learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                       d.phi = 0.001, d.zeta = 0.001,
                       e.phi = 0.001, e.zeta = 0.001,
                       type = "rw-t1", c = 1, L0 = 1
                     ),
                     mcmc.opt = list(
                       chain.length = 12000, burnin = 2000, thin = 10, asis = TRUE
                     ),
                     settings.NegBin = list(
                       alpha.r = 2, beta.r = 1, expansion.steps = 10,
                       width = 1, p.overrelax = 0, accuracy.overrelax = 10
                     ),
                     HPD.coverage = 0.95,
                     R.WAIC = 20,
                     posterior.predictive.matrix = FALSE,
                     random.effects = TRUE,
                     progress.bar = FALSE
){

  prior.reg_nb <- prior.reg_zinb.count
  prior.load_nb <- prior.load_zinb.count
  prior.reg_logit <- prior.reg_zinb.inflation
  prior.load_logit <- prior.load_zinb.inflation

  # model aliases
  if(is.null(model) || length(model) != 1)
    stop("Argument 'model' must be either 'Gaussian', 'Normal', 'Probit', 'Logit', 'NegBin' or 'ZINB'.")
  model <- switch(tolower(model),
                  gaussian = "Gaussian",
                  normal   = "Gaussian",
                  probit   = "Probit",
                  logit    = "Logit",
                  negbin   = "NegBin",
                  zinb     = "ZINB",
                  NULL
  )

  # input checks
  check.panelTVP(formula, data, id, t, model, prior.reg, prior.var, prior.load,
                 prior.reg_nb, prior.load_nb, prior.reg_logit, prior.load_logit,
                 mcmc.opt, settings.NegBin, HPD.coverage, R.WAIC, posterior.predictive.matrix,
                 random.effects, progress.bar)

  # ordering dataset and removing gaps in id
  data$id <- as.numeric(factor(id))
  data$t <- t
  data <- data[order(data$t, data$id),]

  if(model != "ZINB"){

    if(prior.reg$type == "rw-t0") prior.reg$type <- "rw1"
    if(prior.reg$type == "rw-t1") prior.reg$type <- "rw2"
    if(prior.load$type == "rw-t0") prior.load$type <- "rw1"
    if(prior.load$type == "rw-t1") prior.load$type <- "rw2"

    # Gaussian, Probit, Logit and Negative Binomial model ------------------------

    # fit the model
    result <- fit_panelTVP(formula = formula,
                           data = data,
                           model = model,
                           prior.reg = prior.reg,
                           prior.var = prior.var,
                           prior.load = prior.load,
                           mcmc.opt = mcmc.opt,
                           settings.NegBin = settings.NegBin,
                           HPD.coverage = HPD.coverage,
                           random.effects = random.effects,
                           progress.bar = progress.bar)
    # if not random effects structure requested, delete placeholders
    if(!random.effects){
      result$fmean <- NULL
      result$mcmc <- result$mcmc[, !startsWith(colnames(result$mcmc), "lambda")]
      result$mcmc <- result$mcmc[, !(colnames(result$mcmc) %in% c(
        "psi", "phi2", "zeta2",
        "a.phi", "kappa.phi",
        "a.zeta", "kappa.zeta"))
        ]
      result$posterior <- result$posterior[!startsWith(rownames(result$posterior), "lambda"),]
      result$posterior <- result$posterior[!(rownames(result$posterior) %in%
                                             c("abs(psi)", "phi2", "zeta2",
                                               "a.phi", "kappa.phi",
                                               "a.zeta", "kappa.zeta")),]
    }
    # add WAIC, compute posterior predictive and remove chain of factor scores to save memory
    cat(" Computing WAIC ...")
    result$WAIC <- compute_waic(result, random.effects, R.WAIC)
    if(posterior.predictive.matrix){
      if(random.effects){
        result$posterior.predictive <- compute_fitted_Gaussian_Probit_Logit_NegBin(result)
      } else{
        result$posterior.predictive <- compute_fitted_Gaussian_Probit_Logit_NegBin_no.fac(result)
      }
    }
    result$fmcmc <- NULL

    # adding learning settings to output
    hyperpara <- c("a.xi", "a.tau", "c.xi", "c.tau", "kappa.xi", "kappa.tau", "kappa.zeta", "kappa.phi")
    part <- c(rep("regression part", 6), rep("factor part", 2))
    learn <- c(prior.reg$learn.a.xi, prior.reg$learn.a.tau,
               prior.reg$learn.c.xi, prior.reg$learn.c.tau,
               prior.reg$learn.kappa.xi, prior.reg$learn.kappa.tau,
               prior.load$learn.kappa.zeta, prior.load$learn.kappa.phi)
    if(!(prior.reg$type %in% c("rw1", "rw2"))) learn[1:6] <- NA
    if(!(prior.load$type %in% c("rw1", "rw2"))) learn[7:8] <- NA
    result$learning.settings <- cbind(hyperpara, part, learn)
    colnames(result$learning.settings) <- c("hyperparameter", "model.part", "learned?")
    result$learning.settings <- as.data.frame(result$learning.settings)
    if(!random.effects) result$learning.settings <- result$learning.settings[1:6,]
    if(!prior.reg$TG){
      result$learning.settings[3:4, 3] <- NA
    }
    if(prior.reg$TG && prior.reg$TG.alternative){
      result$learning.settings[1:6, 3] <- NA
    }
    # adding mcmc setting to output (incl. ASIS Boolean)
    result$mcmc.settings <- mcmc.opt
    # adding matrix to match effect ids to variable names (for user information)
    x <- result$data$X
    variable.codes <- matrix(nrow = ncol(x), ncol = 2)
    colnames(variable.codes) <- c("variable", "effect")
    variable.codes[,1] <- colnames(x)
    variable.codes[,2] <- paste0("beta",1:ncol(x))
    result$variable.codes <- variable.codes

  } else{

    # Zero-Inflated Negative Binomial model ------------------------------------

    if(prior.reg_nb$type == "rw-t0") prior.reg_nb$type <- "rw1"
    if(prior.reg_nb$type == "rw-t1") prior.reg_nb$type <- "rw2"
    if(prior.load_nb$type == "rw-t0") prior.load_nb$type <- "rw1"
    if(prior.load_nb$type == "rw-t1") prior.load_nb$type <- "rw2"

    if(prior.reg_logit$type == "rw-t0") prior.reg_logit$type <- "rw1"
    if(prior.reg_logit$type == "rw-t1") prior.reg_logit$type <- "rw2"
    if(prior.load_logit$type == "rw-t0") prior.load_logit$type <- "rw1"
    if(prior.load_logit$type == "rw-t1") prior.load_logit$type <- "rw2"

    result <- fit_panelTVP_ZINB(formula = formula,
                                data = data,
                                prior.reg_nb = prior.reg_nb,
                                prior.load_nb = prior.load_nb,
                                prior.reg_logit = prior.reg_logit,
                                prior.load_logit = prior.load_logit,
                                mcmc.opt = mcmc.opt,
                                settings.NegBin = settings.NegBin,
                                HPD.coverage = HPD.coverage,
                                random.effects = random.effects,
                                progress.bar = progress.bar)
    # if not random effects structure requested, delete placeholders
    if(!random.effects){
      result$fmean_nb <- NULL
      result$fmean_logit <- NULL
      result$mcmc_nb <- result$mcmc_nb[, !startsWith(colnames(result$mcmc_nb), "lambda")]
      result$mcmc_logit <- result$mcmc_logit[, !startsWith(colnames(result$mcmc_logit), "lambda")]
      result$mcmc_nb <- result$mcmc_nb[, !(colnames(result$mcmc_nb) %in% c(
        "psi", "phi2", "zeta2",
        "a.phi", "kappa.phi",
        "a.zeta", "kappa.zeta"))
      ]
      result$mcmc_logit <- result$mcmc_logit[, !(colnames(result$mcmc_logit) %in% c(
        "psi", "phi2", "zeta2",
        "a.phi", "kappa.phi",
        "a.zeta", "kappa.zeta"))
      ]
      result$posterior_nb <- result$posterior_nb[!startsWith(rownames(result$posterior_nb), "lambda"),]
      result$posterior_nb <- result$posterior_nb[!(rownames(result$posterior_nb) %in%
                                               c("abs(psi)", "phi2", "zeta2",
                                                 "a.phi", "kappa.phi",
                                                 "a.zeta", "kappa.zeta")),]
      result$posterior_logit <- result$posterior_logit[!startsWith(rownames(result$posterior_logit), "lambda"),]
      result$posterior_logit <- result$posterior_logit[!(rownames(result$posterior_logit) %in%
                                                     c("abs(psi)", "phi2", "zeta2",
                                                       "a.phi", "kappa.phi",
                                                       "a.zeta", "kappa.zeta")),]
    }
    # add WAIC, compute posterior predictive and remove chain of factor scores and risk-indicators to save memory
    cat(" Computing WAIC ...")
    result$WAIC <- compute_waic(result, random.effects, R.WAIC)
    if(posterior.predictive.matrix){
      if(random.effects){
        result$posterior.predictive <- compute_fitted_ZINB(result)
      } else{
        result$posterior.predictive <- compute_fitted_ZINB_no.fac(result)
      }
    }
    result$fmcmc_logit <- NULL
    result$fmcmc_nb <- NULL
    result$mcmc_risk <- NULL

    # adding learning settings to output
    hyperpara <- c("a.xi", "a.tau", "c.xi", "c.tau", "kappa.xi", "kappa.tau", "kappa.zeta", "kappa.phi")
    part <- c(rep("regression part", 6), rep("factor part", 2))
    learn_nb <- c(prior.reg_nb$learn.a.xi, prior.reg_nb$learn.a.tau,
                  prior.reg_nb$learn.c.xi, prior.reg_nb$learn.c.tau,
                  prior.reg_nb$learn.kappa.xi, prior.reg_nb$learn.kappa.tau,
                  prior.load_nb$learn.kappa.zeta, prior.load_nb$learn.kappa.phi)
    if(!(prior.reg_nb$type %in% c("rw1", "rw2"))) learn_nb[1:6] <- NA
    if(!(prior.load_nb$type %in% c("rw1", "rw2"))) learn_nb[7:8] <- NA
    learn_logit <- c(prior.reg_logit$learn.a.xi, prior.reg_logit$learn.a.tau,
                     prior.reg_logit$learn.c.xi, prior.reg_logit$learn.c.tau,
                     prior.reg_logit$learn.kappa.xi, prior.reg_logit$learn.kappa.tau,
                     prior.load_logit$learn.kappa.zeta, prior.load_logit$learn.kappa.phi)
    if(!(prior.reg_logit$type %in% c("rw1", "rw2"))) learn_logit[1:6] <- NA
    if(!(prior.load_logit$type %in% c("rw1", "rw2"))) learn_logit[7:8] <- NA
    result$learning.settings_logit <- cbind(hyperpara, part, learn_logit)
    colnames(result$learning.settings_logit) <- c("hyperparameter", "model.part", "learned?")
    result$learning.settings_logit <- as.data.frame(result$learning.settings_logit)
    result$learning.settings_nb <- cbind(hyperpara, part, learn_nb)
    colnames(result$learning.settings_nb) <- c("hyperparameter", "model.part", "learned?")
    result$learning.settings_nb <- as.data.frame(result$learning.settings_nb)
    if(!random.effects){
      result$learning.settings_nb <- result$learning.settings_nb[1:6,]
      result$learning.settings_logit <- result$learning.settings_logit[1:6,]
    }
    if(!prior.reg_nb$TG){
      result$learning.settings_nb[3:4, 3] <- NA
    }
    if(!prior.reg_logit$TG){
      result$learning.settings_nb[3:4, 3] <- NA
    }
    if(prior.reg_nb$TG && prior.reg_nb$TG.alternative){
      result$learning.settings_nb[1:6, 3] <- NA
    }
    if(prior.reg_logit$TG && prior.reg_logit$TG.alternative){
      result$learning.settings_logit[1:6, 3] <- NA
    }
    # adding mcmc setting to output (incl. ASIS Boolean)
    result$mcmc.settings <- mcmc.opt

    # adding matrix to match effect ids to variable names (for user information)
    x_nb <- result$data$X_nb
    variable.codes_nb <- matrix(nrow = ncol(x_nb), ncol = 2)
    colnames(variable.codes_nb) <- c("variable", "effect")
    variable.codes_nb[,1] <- colnames(x_nb)
    variable.codes_nb[,2] <- paste0("beta",1:ncol(x_nb))
    result$variable.codes_nb <- as.data.frame(variable.codes_nb)

    x_logit <- result$data$X_logit
    variable.codes_logit <- matrix(nrow = ncol(x_logit), ncol = 2)
    colnames(variable.codes_logit) <- c("variable", "effect")
    variable.codes_logit[,1] <- colnames(x_logit)
    variable.codes_logit[,2] <- paste0("beta",1:ncol(x_logit))
    result$variable.codes_logit <- as.data.frame(variable.codes_logit)

  }

  cat(" Analysis finished!")
  return(result)

}
