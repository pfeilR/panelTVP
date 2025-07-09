#' @title Print basic model information output
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.Gaussian}.
#'
#' @param x an object of class \code{panelTVP.Gaussian}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @exportS3Method print panelTVP.Gaussian
#' @examples
#' # Printing object of class panelTVP.Gaussian
#' # NB: To reduces computational effort, we have drastically reduced the length
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
#'                          mcmc.opt = list(chain.length = 500, burnin = 100, thin = 1, asis = TRUE),
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
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.Probit}.
#'
#' @param x an object of class \code{panelTVP.Probit}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @exportS3Method print panelTVP.Probit
#'
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
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.Logit}.
#'
#' @param x an object of class \code{panelTVP.Logit}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @exportS3Method print panelTVP.Logit
#'
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
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.NegBin}.
#'
#' @param x an object of class \code{panelTVP.NegBin}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @exportS3Method print panelTVP.NegBin
#'
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
  - learning.settings: Information on learning status of hyperparameters. \n
  - mcmc.settings: Details on MCMC sampler. \n
You may use the following functions to get additional information:\n
  - summary(): This will give you a formatted summary of the most important parameters.\n
  - plot(): This will give you plots of the coefficient estimates.\n
  - predict(): This will give you the posterior predictive distribution of new data.\n
")
  invisible(x)
}

#' @title Print basic model information output
#'
#' @description
#'  This basic \code{print} method gives a general overview
#'    of the information that is contained in an object of class \code{panelTVP.ZINB}.
#'
#' @param x an object of class \code{panelTVP.ZINB}
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @exportS3Method print panelTVP.ZINB
#'
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
