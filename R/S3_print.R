#' @title Print basic model information output for a \code{panelTVP.Gaussian} object
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.Gaussian}.
#'
#' @param x an object of class \code{panelTVP.Gaussian}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @name print.panelTVP.Gaussian
#' @rdname print.panelTVP.Gaussian
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method print panelTVP.Gaussian
#' @examples
#' # Printing object of class panelTVP.Gaussian
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.gaussian <- sim_panelTVP(n = 100,
#'                              Tmax = 4,
#'                              beta = c(4,1,0,0),
#'                              theta = c(1,0.5,0,0),
#'                              lambda = 1,
#'                              psi = 0.2,
#'                              model = "Gaussian",
#'                              sigma2 = 0.7)
#' res.gaussian <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#' print(res.gaussian)
print.panelTVP.Gaussian <- function(x, ...){
  cat("\nThis is an object of class panelTVP.Gaussian. It contains:\n
  - data: Input data and additional context information derived from the data.\n
  - mcmc: Markov Chains for every parameter except for the factor scores.\n
  - posterior: Posterior summary. \n
  - fmean: Posterior means of random intercepts.\n
  - model: The model you have fitted.\n
  - acceptance.rates: The achieved acceptance rates of Metropolis-Hastings.\n
  - HPD.coverage: Coverage probability of Highest Posterior Density Intervals. \n
  - runtime: The total time for fitting the model (measured in seconds).\n
  - WAIC: The Widely Applicable Information Criterion (or Watanabe's AIC).\n
  - fitted.values: The fitted values for each observation (rows) and each MCMC draw (columns). \n
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output for a \code{panelTVP.Probit} object
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.Probit}.
#'
#' @param x an object of class \code{panelTVP.Probit}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method print panelTVP.Probit
#' @examples
#' # Printing object of class panelTVP.Probit
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.probit <- sim_panelTVP(n = 100,
#'                            Tmax = 4,
#'                            beta = c(1,0.5,0,0),
#'                            theta = c(0.8,0.5,0,0),
#'                            lambda = 1,
#'                            psi = 0.2,
#'                            model = "Probit")
#' res.probit <- panelTVP(y ~ W1 + W2 + W3,
#'                        data = sim.probit$observed,
#'                        mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                        model = "Probit")
#' print(res.probit)
print.panelTVP.Probit <- function(x, ...){
  cat("\nThis is an object of class panelTVP.Probit. It contains:\n
  - data: Input data and additional context information derived from the data.\n
  - mcmc: Markov Chains for every parameter except for the factor scores.\n
  - posterior: Posterior summary. \n
  - fmean: Posterior means of random intercepts.\n
  - model: The model you have fitted.\n
  - acceptance.rates: The achieved acceptance rates of Metropolis-Hastings.\n
  - HPD.coverage: Coverage probability of Highest Posterior Density Intervals. \n
  - runtime: The total time for fitting the model (measured in seconds).\n
  - WAIC: The Widely Applicable Information Criterion (or Watanabe's AIC).\n
  - fitted.values: The fitted values for each observation (rows) and each MCMC draw (columns). \n
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output for a \code{panelTVP.Logit} object
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.Logit}.
#'
#' @param x an object of class \code{panelTVP.Logit}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method print panelTVP.Logit
#' @examples
#' # Printing object of class panelTVP.Logit
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.logit <- sim_panelTVP(n = 100,
#'                           Tmax = 4,
#'                           beta = c(1,0.5,0,0),
#'                           theta = c(0.8,0.5,0,0),
#'                           lambda = 1,
#'                           psi = 0.2,
#'                           model = "Logit")
#' res.logit <- panelTVP(y ~ W1 + W2 + W3,
#'                       data = sim.logit$observed,
#'                       mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                       model = "Logit")
#' print(res.logit)
print.panelTVP.Logit <- function(x, ...){
  cat("\nThis is an object of class panelTVP.Logit. It contains:\n
  - data: Input data and additional context information derived from the data.\n
  - mcmc: Markov Chains for every parameter except for the factor scores.\n
  - posterior: Posterior summary. \n
  - fmean: Posterior means of random intercepts.\n
  - model: The model you have fitted.\n
  - acceptance.rates: The achieved acceptance rates of Metropolis-Hastings.\n
  - HPD.coverage: Coverage probability of Highest Posterior Density Intervals. \n
  - runtime: The total time for fitting the model (measured in seconds).\n
  - WAIC: The Widely Applicable Information Criterion (or Watanabe's AIC).\n
  - fitted.values: The fitted values for each observation (rows) and each MCMC draw (columns). \n
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output for a \code{panelTVP.NegBin} object
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.NegBin}.
#'
#' @param x an object of class \code{panelTVP.NegBin}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method print panelTVP.NegBin
#' @examples
#' # Printing object of class panelTVP.NegBin
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.negbin <- sim_panelTVP(n = 100,
#'                            Tmax = 4,
#'                            beta = c(1,0.5,0,0),
#'                            theta = c(0.8,0.5,0,0),
#'                            lambda = 1,
#'                            psi = 0.2,
#'                            r = 2,
#'                            model = "NegBin")
#' res.negbin <- panelTVP(y ~ W1 + W2 + W3,
#'                        data = sim.negbin$observed,
#'                        mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                        model = "NegBin")
#' print(res.negbin)
print.panelTVP.NegBin <- function(x, ...){
  cat("\nThis is an object of class panelTVP.NegBin. It contains:\n
  - data: Input data and additional context information derived from the data.\n
  - mcmc: Markov Chains for every parameter except for the factor scores.\n
  - posterior: Posterior summary. \n
  - fmean: Posterior means of random intercepts.\n
  - model: The model you have fitted.\n
  - acceptance.rates: The achieved acceptance rates of Metropolis-Hastings.\n
  - HPD.coverage: Coverage probability of Highest Posterior Density Intervals. \n
  - runtime: The total time for fitting the model (measured in seconds).\n
  - WAIC: The Widely Applicable Information Criterion (or Watanabe's AIC).\n
  - fitted.values: The fitted values for each observation (rows) and each MCMC draw (columns). \n
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output for a \code{panelTVP.ZINB} object
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.ZINB}.
#'
#' @param x an object of class \code{panelTVP.ZINB}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method print panelTVP.ZINB
#' @examples
#' # Printing object of class panelTVP.ZINB
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.zinb <- sim_panelTVP(n = 100,
#'                          Tmax = 4,
#'                          beta.nb = c(0.5,-0.7,0,0),
#'                          theta.nb = c(0.05,0.5,0,0),
#'                          lambda.nb = 0.5,
#'                          psi.nb = 0.02,
#'                          beta.logit = c(-1,0.6,0,0),
#'                          theta.logit = c(0,1,0,0),
#'                          lambda.logit = 0.7,
#'                          psi.logit = 0,
#'                          r = 2,
#'                          model = "ZINB")
#' res.zinb <- panelTVP(y ~ W1.nb + W2.nb + W3.nb | W1.logit + W2.logit + W3.logit,
#'                      data = sim.zinb$observed,
#'                      mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                      model = "ZINB")
#' print(res.zinb)
print.panelTVP.ZINB <- function(x, ...){
  cat("\nThis is an object of class panelTVP.ZINB. It contains:\n
  - data: Input data and additional context information derived from the data.\n
  - mcmc_logit: Markov Chains for every parameter except for the factor scores. (Logit) \n
  - mcmc_nb: Markov Chains for every parameter except for the factor scores. (Negative Binomial) \n
  - posterior_logit: Posterior summary. (Logit) \n
  - posterior_nb: Posterior summary. (Negative Binomial) \n
  - fmean_logit: Posterior means of random intercepts. (Logit) \n
  - fmean_nb: Posterior means of random intercepts. (Negative Binomial) \n
  - model: The model you have fitted.\n
  - acceptance.rates: The achieved acceptance rates of Metropolis-Hastings.\n
  - HPD.coverage: Coverage probability of Highest Posterior Density Intervals. \n
  - runtime: The total time for fitting the model (measured in seconds).\n
  - WAIC: The Widely Applicable Information Criterion (or Watanabe's AIC).\n
  - fitted.values: The fitted values for each observation (rows) and each MCMC draw (columns). \n
  - learning.settings_logit: Information on learning status of hyperparameters. (Logit) \n
  - learning.settings_nb: Information on learning status of hyperparameters. (Negative Binomial) \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}
