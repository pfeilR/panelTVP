#' Fit a Bayesian IV panel data model with time-varying parameters
#'
#' @usage panelTVP_IV(formula_stage1 = NULL,
#'             formula_stage2 = NULL,
#'             data = NULL,
#'             id = NULL,
#'             t = NULL,
#'             prior.reg_stage1 = list(),
#'             prior.reg_stage2 = list(),
#'             prior.var_stage2 = list(),
#'             prior.load_stage2 = list(),
#'             prior.rho = list(),
#'             mcmc.opt = list(),
#'             HPD.coverage = 0.95,
#'             posterior.predictive.matrix = FALSE,
#'             random.effects = TRUE,
#'             progress.bar = FALSE)
#'
#' @param formula_stage1 formula object for the first-stage Probit model (no default)
#' @param formula_stage2 formula object for the second-stage Gaussian model (no default)
#' @param data a data frame that contains the variables of the formula arguments
#'  (no default)
#' @param id a vector with length equal to the number of observations in the dataset, i.e.,
#'  this is the 'subject' variable in the dataset (no default)
#' @param t a vector with length equal to the number of observations in the dataset, i.e.,
#'  this is the 'time' variable in the dataset (no default)
#' @param prior.reg_stage1 a list of arguments for estimating the parameters of the regression
#'  part of the first-stage Probit model. Note that this stage does not contain
#'  random effects. The prior arguments are the ones of [panelTVP()]
#' @param prior.reg_stage2 a list of arguments for estimating the parameters of the regression
#'  part of the second-stage Gaussian model. The prior arguments are the ones of [panelTVP()].
#' @param prior.var_stage2 a list of arguments for estimating the error variance
#'  of the second-stage Gaussian model using 2D-Slice sampling. The arguments are:
#'   \itemize{
#'    \item \code{alpha.sigma}: shape parameter of Inverse-Gamma prior
#'    \item \code{beta.sigma}: rate parameter of Inverse-Gamma prior
#'    \item \code{expansion.steps}: number of steps in stepping-out phase
#'    \item \code{width}: width of the slice interval (on log-scale of error variance)
#'   }
#' @param prior.load_stage2 a list of arguments for estimating the parameters of the factor
#'  part of the second-stage Gaussian model. The prior arguments are the ones of [panelTVP()]
#' @param prior.rho a list of arguments for estimating the correlation of errors,
#'  i.e. the degree of endogeneity using 2D-Slice sampling. The arguments are:
#'   \itemize{
#'    \item \code{mean.rho}: mean of Normal prior on Fisher-Z scale of \code{rho}
#'    \item \code{sd.rho}: standard deviation of Normal prior on Fisher-Z scale of \code{rho}
#'    \item \code{expansion.steps}: number of steps in stepping-out phase
#'    \item \code{width}: width of the slice interval (on Fisher-Z scale of \code{rho})
#'   }
#' @param mcmc.opt a list containing information on the overall sampler.
#'  The arguments are:
#'  \itemize{
#'   \item \code{chain.length}: the length of the Markov Chain (i.e., total number
#'   of draws)
#'   \item \code{burnin}: the burn-in period
#'   \item \code{thin}: the thinning factor
#'   \item \code{asis}: if set to TRUE, an ancillarity sufficiency interweaving
#'    step is added for increasing the sampling efficiency of the regression
#'    effects and factor loadings
#'  }
#' @param HPD.coverage coverage probability of highest posterior density intervals
#'  (default yields 95 percent coverage)
#' @param posterior.predictive.matrix if TRUE (= default) the posterior predictive distribution
#'  based on the training data is computed and returned as a list object by the function.
#'  Setting it to FALSE is usually only done for saving memory
#' @param random.effects if TRUE (= default) a factor model is included for estimating
#'  random effects in the second-stage of the IV model,
#'  if FALSE the model does not contain random effects and, consequently,
#'  priors on the parameters of the factor model are ignored
#' @param progress.bar if TRUE a progress bar is displayed, if FALSE the progress bar is omitted
#'
#' @description
#' This function estimates a Bayesian Instrumental Variable Model for identifying
#'  causal treatment effects. The function assumes that the (endogenous) treatment
#'  is binary whereas a Gaussian model is used for the observation equation (second stage).
#'  Similarly to the function [panelTVP()], effect estimates vary over time and
#'  shrinkage priors are used to classify effects to time-varying, constant and
#'  zero. Unobserved heterogeneity is captured by random effects at the second stage
#'  via a factor model. The error variance of the second stage model as well as
#'  the correlation between errors of both stages are sampled jointly via
#'  2D-Slice sampling (Neal, 2003). For this, an Inverse-Gamma prior is used for the variance
#'  whereas a Normal prior is placed on the Fisher-Z scale of the correlation
#'  coefficient. Time-varying effects are modelled with the same priors as in
#'  the non-treatment-effects models in [panelTVP()].
#'
#' @returns The function returns an object of class \code{panelTVP.IV} and
#' contains a list of the following elements:
#'  \describe{
#'    \item{learning.settings_stage1}{information on which parameters have been learned
#'       for the first stage of the model}
#'    \item{learning.settings_stage2}{information on which parameters have been learned
#'       for the second stage of the model}
#'    \item{variable.codes_stage1}{information on which covariate is associated with which parameter
#'       for the first stage of the model}
#'    \item{variable.codes_stage2}{information on which covariate is associated with which parameter
#'       for the second stage of the model}
#'    \item{data}{the data used for fitting the model and additional context information
#'    derived from the data}
#'    \item{Y}{the \eqn{Tn \times M} response data matrix of every iteration of the chain,
#'    i.e., this matrix is only included in the output when missing response data
#'    were present and imputed via data augmentation. Each row of \code{Y} contains
#'    the sampled values of a specific observation. For observed data, the corresponding
#'    rows are essentially replicates of the same value. For missing data, the
#'    corresponding rows contain the imputed values of every iteration.}
#'    \item{mcmc_stage1}{Markov Chains for every parameter except for the factor scores
#'     (to save memory) for the first stage of the model}
#'    \item{mcmc_stage2}{Markov Chains for every parameter except for the factor scores
#'     (to save memory) for the second stage of the model}
#'    \item{mcmc_rho}{a vector of dimension M containing the Markov Chain of the
#'     correlation coefficient}
#'    \item{posterior_stage1}{preliminary summary of posterior results for the
#'     first stage of the model}
#'    \item{posterior_stage2}{preliminary summary of posterior results for the
#'     second stage of the model}
#'    \item{posterior_rho}{preliminary summary of posterior results for the
#'     correlation coefficient}
#'    \item{fmean_stage2}{posterior means of factor scores for the
#'      second stage of the model (only when random effects are estimated)}
#'    \item{model}{the fitted model}
#'    \item{acceptance.rates}{the achieved acceptance rates when using
#'       Metropolis-Hastings for both components of the model}
#'    \item{HPD.coverage}{coverage probability of HPD intervals (based on input)}
#'    \item{runtime}{total runtime of the sampler (measured in seconds)}
#'    \item{Treatment.Variable}{name of the treatment variable of the model}
#'    \item{posterior.predictive}{the \eqn{Tn \times M} matrix that contains the posterior
#'     predictive distribution for each observation (rows) and each MCMC draw (columns)}
#'    \item{mcmc.settings}{details on the MCMC sampler}
#'  }
#' @references
#'  Neal, R.M. (2003). Slice sampling. In: The Annals of Statistics, 31, 705-767.
#' @examples
#' # Fitting Bayesian IV panel model with time-varying parameters
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.iv <- sim_panelTVP_IV(n = 1000,
#'                           Tmax = 4,
#'                           beta_stage1 = c(1, 0.7),
#'                           theta_stage1 = c(0.2, 0.01),
#'                           beta_stage2 = c(0, 4, 1),
#'                           theta_stage2 = c(0, 3, 0),
#'                           lambda_stage2 = 1.3,
#'                           psi_stage2 = 0.1,
#'                           beta_D = 2,
#'                           theta_D = 0.7,
#'                           rho = 0.1,
#'                           sigma2 = 1,
#'                           n.instruments = 1)
#' res.iv <- panelTVP_IV(formula_stage1 = D ~ X_stage1.Z1,
#'                       formula_stage2 = y ~ X_stage2.W1 + X_stage2.W2 + D,
#'                       data = sim.iv$observed,
#'                       id = sim.iv$observed$id,
#'                       t = sim.iv$observed$t,
#'                       prior.rho = list(
#'                        alpha.rho = 1, beta.rho = 1,
#'                        expansion.steps = 10, width = 0.1
#'                       ),
#'                       mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE))
#' @export
panelTVP_IV <- function(formula_stage1 = NULL,
                        formula_stage2 = NULL,
                        data = NULL,
                        id = NULL,
                        t = NULL,
                        prior.reg_stage1 = list(
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
                        prior.reg_stage2 = list(
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
                        prior.var_stage2 = list(
                          alpha.sigma = 2, beta.sigma = 1,
                          expansion.steps = 10, width = 0.1
                        ),
                        prior.load_stage2 = list(
                          a.phi = 0.1, a.zeta = 0.1,
                          kappa.phi = 10, kappa.zeta = 10,
                          learn.kappa.phi = TRUE, learn.kappa.zeta = TRUE,
                          d.phi = 0.001, d.zeta = 0.001,
                          e.phi = 0.001, e.zeta = 0.001,
                          type = "rw-t1", c = 1, L0 = 1
                        ),
                        prior.rho = list(
                          alpha.rho = 1, beta.rho = 1,
                          expansion.steps = 10, width = 0.1
                        ),
                        mcmc.opt = list(
                          chain.length = 12000, burnin = 2000, thin = 10, asis = TRUE
                          ),
                        HPD.coverage = 0.95,
                        posterior.predictive.matrix = FALSE,
                        random.effects = TRUE,
                        progress.bar = FALSE
){

  # HERE INPUT CHECKS ARE NEEDED !!!

  # ordering dataset and removing gaps in id
  data$id <- as.numeric(factor(id))
  data$t <- t
  data <- data[order(data$t, data$id),]

  if(prior.reg_stage1$type == "rw-t0") prior.reg_stage1$type <- "rw1"
  if(prior.reg_stage1$type == "rw-t1") prior.reg_stage1$type <- "rw2"
  if(prior.reg_stage2$type == "rw-t0") prior.reg_stage2$type <- "rw1"
  if(prior.reg_stage2$type == "rw-t1") prior.reg_stage2$type <- "rw2"
  if(prior.load_stage2$type == "rw-t0") prior.load_stage2$type <- "rw1"
  if(prior.load_stage2$type == "rw-t1") prior.load_stage2$type <- "rw2"

  result <- fit_panelTVP_IV(formula_stage1 = formula_stage1,
                            formula_stage2 = formula_stage2,
                            data = data,
                            prior.reg_stage1 = prior.reg_stage1,
                            prior.reg_stage2 = prior.reg_stage2,
                            prior.load_stage2 = prior.load_stage2,
                            prior.var_stage2 = prior.var_stage2,
                            prior.rho = prior.rho,
                            mcmc.opt = mcmc.opt,
                            HPD.coverage = HPD.coverage,
                            random.effects = random.effects,
                            progress.bar = progress.bar)

  # if not random effects structure requested, delete placeholders
  if(!random.effects){
    result$fmean_stage2 <- NULL
    result$mcmc_stage2 <- result$mcmc_stage2[, !startsWith(colnames(result$mcmc_stage2), "lambda")]
    result$mcmc_stage2 <- result$mcmc_stage2[, !(colnames(result$mcmc_stage2) %in% c(
      "psi", "phi2", "zeta2",
      "a.phi", "kappa.phi",
      "a.zeta", "kappa.zeta"))
    ]
    result$posterior_stage2 <- result$posterior_stage2[!startsWith(rownames(result$posterior_stage2), "lambda"),]
    result$posterior_stage2 <- result$posterior_stage2[!(rownames(result$posterior_stage2) %in%
                                                   c("abs(psi)", "phi2", "zeta2",
                                                     "a.phi", "kappa.phi",
                                                     "a.zeta", "kappa.zeta")),]
  }

  # posterior predictive distribution
  if(posterior.predictive.matrix){
    if(random.effects){
      result$posterior.predictive <- compute_fitted_IV(result)
    } else{
      result$posterior.predictive <- compute_fitted_IV_no.fac(result)
    }
  }
  result$D_mcmc <- NULL
  result$fmcmc_stage2 <- NULL

  # adding learning settings to output
  hyperpara <- c("a.xi", "a.tau", "c.xi", "c.tau", "kappa.xi", "kappa.tau", "kappa.zeta", "kappa.phi")
  part <- c(rep("regression part", 6), rep("factor part", 2))
  learn_stage2 <- c(prior.reg_stage2$learn.a.xi, prior.reg_stage2$learn.a.tau,
                prior.reg_stage2$learn.c.xi, prior.reg_stage2$learn.c.tau,
                prior.reg_stage2$learn.kappa.xi, prior.reg_stage2$learn.kappa.tau,
                prior.load_stage2$learn.kappa.zeta, prior.load_stage2$learn.kappa.phi)
  if(!(prior.reg_stage2$type %in% c("rw1", "rw2"))) learn_stage2[1:6] <- NA
  if(!(prior.load_stage2$type %in% c("rw1", "rw2"))) learn_stage2[7:8] <- NA
  learn_stage1 <- c(prior.reg_stage1$learn.a.xi, prior.reg_stage1$learn.a.tau,
                   prior.reg_stage1$learn.c.xi, prior.reg_stage1$learn.c.tau,
                   prior.reg_stage1$learn.kappa.xi, prior.reg_stage1$learn.kappa.tau)
  if(!(prior.reg_stage1$type %in% c("rw1", "rw2"))) learn_stage1[1:6] <- NA
  result$learning.settings_stage1 <- cbind(hyperpara[1:6], part[1:6], learn_stage1)
  colnames(result$learning.settings_stage1) <- c("hyperparameter", "model.part", "learned?")
  result$learning.settings_stage1 <- as.data.frame(result$learning.settings_stage1)
  result$learning.settings_stage2 <- cbind(hyperpara, part, learn_stage2)
  colnames(result$learning.settings_stage2) <- c("hyperparameter", "model.part", "learned?")
  result$learning.settings_stage2 <- as.data.frame(result$learning.settings_stage2)
  if(!random.effects){
    result$learning.settings_stage1 <- result$learning.settings_stage1[1:6,]
    result$learning.settings_stage2 <- result$learning.settings_stage2[1:6,]
  }
  if(!prior.reg_stage1$TG){
    result$learning.settings_stage1[3:4, 3] <- NA
  }
  if(!prior.reg_stage2$TG){
    result$learning.settings_stage2[3:4, 3] <- NA
  }
  # adding mcmc setting to output (incl. ASIS Boolean)
  result$mcmc.settings <- mcmc.opt

  # adding matrix to match effect ids to variable names (for user information)
  x_stage1 <- result$data$X_stage1
  variable.codes_stage1 <- matrix(nrow = ncol(x_stage1), ncol = 2)
  colnames(variable.codes_stage1) <- c("variable", "effect")
  variable.codes_stage1[,1] <- colnames(x_stage1)
  variable.codes_stage1[,2] <- paste0("beta",1:ncol(x_stage1))
  result$variable.codes_stage1 <- as.data.frame(variable.codes_stage1)

  x_stage2 <- result$data$X_stage2
  variable.codes_stage2 <- matrix(nrow = ncol(x_stage2), ncol = 2)
  colnames(variable.codes_stage2) <- c("variable", "effect")
  variable.codes_stage2[,1] <- colnames(x_stage2)
  variable.codes_stage2[,2] <- paste0("beta",1:ncol(x_stage2))
  result$variable.codes_stage2 <- as.data.frame(variable.codes_stage2)

  cat(" Analysis finished!")
  return(result)

}
