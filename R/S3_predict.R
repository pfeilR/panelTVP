#' @title Get predictions for new subjects based on a \code{panelTVP.Gaussian} object
#'
#' @description
#'  This \code{predict} function simulates data and computes summary statistics based on the
#'   posterior predictive distribution for new data under a Gaussian likelihood.
#'
#' @param object an object of class \code{panelTVP.Gaussian}
#' @param X.new a matrix or data frame consisting of new data for the same variables
#'  that were used for fitting the model. The first column must contain a 1, when
#'  the fitted model contains an intercept.
#' @param timepoint a numeric scalar indicating the time point for which predictions should
#'  be made, i.e., predictions for a given data set are only made for one specific time point.
#'  In case you want predictions for additional time points you need to repeatedly call this
#'  function.
#' @param coverage coverage probability for prediction intervals - defaults to 95 percent coverage
#' @param pop.pred logical value, if TRUE population-based predictions are made,
#'  that ignore the random effects structure, if FALSE the random effects structure
#'  is included as well where the unknown factor scores are sampled from their
#'  standard Normal prior - defaults to FALSE
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @returns a list containing the following elements
#' \itemize{
#'   \item \code{predictive.distribution}: a matrix where each row contains draws
#'    from the posterior predictive distribution for the corresponding observation
#'   \item \code{predictive.summary}: posterior mean as well as HPD interval based on
#'    posterior predictive distribution for each new observation
#' }
#'
#' @exportS3Method predict panelTVP.Gaussian
#' @examples
#' # Predictions based on an object of class panelTVP.Gaussian
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
#' # setting up design matrix for predicting two new observations
#' X.new <- data.frame(cbind(c(1,1), # column of 1's for the intercept
#'                           c(2,0),
#'                           c(4,3),
#'                           c(0,0)))
#' colnames(X.new) <- colnames(res.gaussian$data$X)
#' # prediction for the 2nd panel wave
#' pp <- predict(res.gaussian, X.new = X.new, timepoint = 2)
#' plot(density(pp$predictive.distribution[1,]))
#' plot(density(pp$predictive.distribution[2,]))
#' pp$predictive.summary
predict.panelTVP.Gaussian <- function(object, X.new, timepoint,
                                      coverage = 0.95, pop.pred = FALSE, ...){
  check.predict(model = object, X.new = X.new, timepoint = timepoint, coverage = coverage,
                pop.pred = pop.pred)
  pred.helper(model = object, X.new = X.new, timepoint = timepoint,
                 coverage = coverage, pop.pred = pop.pred)
}

#' @title Get predictions for new subjects based on a \code{panelTVP.Probit} object
#'
#' @description
#'  This \code{predict} function simulates data and computes summary statistics based on the
#'   posterior predictive distribution for new data under a Bernoulli likelihood using a
#'   Probit link.
#'
#' @param object an object of class \code{panelTVP.Probit}
#' @param X.new a matrix or data frame consisting of new data for the same variables
#'  that were used for fitting the model. The first column must contain a 1, when
#'  the fitted model contains an intercept.
#' @param timepoint a numeric scalar indicating the time point for which predictions should
#'  be made, i.e., predictions for a given data set are only made for one specific time point.
#'  In case you want predictions for additional time points you need to repeatedly call this
#'  function.
#' @param coverage coverage probability for prediction intervals - defaults to 95 percent coverage
#' @param pop.pred logical value, if TRUE population-based predictions are made,
#'  that ignore the random effects structure, if FALSE the random effects structure
#'  is included as well where the unknown factor scores are sampled from their
#'  standard Normal prior - defaults to FALSE
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @returns a list containing the following elements
#' \itemize{
#'   \item \code{predictive.distribution}: a matrix where each row contains draws
#'    from the posterior predictive distribution for the corresponding observation
#'   \item \code{predictive.summary}: posterior mean as well as HPD interval based on
#'    posterior predictive distribution for each new observation
#' }
#' @exportS3Method predict panelTVP.Probit
#' @examples
#' # Predictions based on an object of class panelTVP.Probit
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
#' # setting up design matrix for predicting two new observations
#' X.new <- data.frame(cbind(c(1,1), # column of 1's for the intercept
#'                           c(2,0),
#'                           c(4,3),
#'                           c(0,0)))
#' colnames(X.new) <- colnames(res.probit$data$X)
#' # prediction for the 2nd panel wave
#' pp <- predict(res.probit, X.new = X.new, timepoint = 2)
#' plot(density(pp$predictive.distribution[1,]))
#' plot(density(pp$predictive.distribution[2,]))
#' pp$predictive.summary
predict.panelTVP.Probit <- function(object, X.new, timepoint,
                                    coverage = 0.95, pop.pred = FALSE, ...){
  check.predict(model = object, X.new = X.new, timepoint = timepoint, coverage = coverage,
                pop.pred = pop.pred)
  pred.helper(model = object, X.new = X.new, timepoint = timepoint,
                 coverage = coverage, pop.pred = pop.pred)
}

#' @title Get predictions for new subjects based on a \code{panelTVP.Logit} object
#'
#' @description
#'  This \code{predict} function simulates data and computes summary statistics based on the
#'   posterior predictive distribution for new data under a Bernoulli likelihood using a
#'   Logit link.
#'
#' @param object an object of class \code{panelTVP.Logit}
#' @param X.new a matrix or data frame consisting of new data for the same variables
#'  that were used for fitting the model. The first column must contain a 1, when
#'  the fitted model contains an intercept.
#' @param timepoint a numeric scalar indicating the time point for which predictions should
#'  be made, i.e., predictions for a given data set are only made for one specific time point.
#'  In case you want predictions for additional time points you need to repeatedly call this
#'  function.
#' @param coverage coverage probability for prediction intervals - defaults to 95 percent coverage
#' @param pop.pred logical value, if TRUE population-based predictions are made,
#'  that ignore the random effects structure, if FALSE the random effects structure
#'  is included as well where the unknown factor scores are sampled from their
#'  standard Normal prior - defaults to FALSE
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @returns a list containing the following elements
#' \itemize{
#'   \item \code{predictive.distribution}: a matrix where each row contains draws
#'    from the posterior predictive distribution for the corresponding observation
#'   \item \code{predictive.summary}: posterior mean as well as HPD interval based on
#'    posterior predictive distribution for each new observation
#' }
#' @exportS3Method predict panelTVP.Logit
#' @examples
#' # Predictions based on an object of class panelTVP.Logit
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
#' # setting up design matrix for predicting two new observations
#' X.new <- data.frame(cbind(c(1,1), # column of 1's for the intercept
#'                           c(2,0),
#'                           c(4,3),
#'                           c(0,0)))
#' colnames(X.new) <- colnames(res.logit$data$X)
#' # prediction for the 2nd panel wave
#' pp <- predict(res.logit, X.new = X.new, timepoint = 2)
#' plot(density(pp$predictive.distribution[1,]))
#' plot(density(pp$predictive.distribution[2,]))
#' pp$predictive.summary
predict.panelTVP.Logit <- function(object, X.new, timepoint,
                                   coverage = 0.95, pop.pred = FALSE, ...){
  check.predict(model = object, X.new = X.new, timepoint = timepoint, coverage = coverage,
                pop.pred = pop.pred)
  pred.helper(model = object, X.new = X.new, timepoint = timepoint,
                 coverage = coverage, pop.pred = pop.pred)
}

#' @title Get predictions for new subjects based on a \code{panelTVP.NegBin} object
#'
#' @description
#'  This \code{predict} function simulates data and computes summary statistics based on the
#'   posterior predictive distribution for new data under a Negative Binomial likelihood.
#'
#' @param object an object of class \code{panelTVP.NegBin}
#' @param X.new a matrix or data frame consisting of new data for the same variables
#'  that were used for fitting the model. The first column must contain a 1, when
#'  the fitted model contains an intercept.
#' @param timepoint a numeric scalar indicating the time point for which predictions should
#'  be made, i.e., predictions for a given data set are only made for one specific time point.
#'  In case you want predictions for additional time points you need to repeatedly call this
#'  function.
#' @param coverage coverage probability for prediction intervals - defaults to 95 percent coverage
#' @param pop.pred logical value, if TRUE population-based predictions are made,
#'  that ignore the random effects structure, if FALSE the random effects structure
#'  is included as well where the unknown factor scores are sampled from their
#'  standard Normal prior - defaults to FALSE
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @returns a list containing the following elements
#' \itemize{
#'   \item \code{predictive.distribution}: a matrix where each row contains draws
#'    from the posterior predictive distribution for the corresponding observation
#'   \item \code{predictive.summary}: posterior mean as well as HPD interval based on
#'    posterior predictive distribution for each new observation
#' }
#' @exportS3Method predict panelTVP.NegBin
#' @examples
#' # Predictions based on an object of class panelTVP.NegBin
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
#' # setting up design matrix for predicting two new observations
#' X.new <- data.frame(cbind(c(1,1), # column of 1's for the intercept
#'                           c(2,0),
#'                           c(4,3),
#'                           c(0,0)))
#' colnames(X.new) <- colnames(res.negbin$data$X)
#' # prediction for the 2nd panel wave
#' pp <- predict(res.negbin, X.new = X.new, timepoint = 2)
#' plot(density(pp$predictive.distribution[1,]))
#' plot(density(pp$predictive.distribution[2,]))
#' pp$predictive.summary
predict.panelTVP.NegBin <- function(object, X.new, timepoint,
                                    coverage = 0.95, pop.pred = FALSE, ...){
  check.predict(model = object, X.new = X.new, timepoint = timepoint, coverage = coverage,
                pop.pred = pop.pred)
  pred.helper(model = object, X.new = X.new, timepoint = timepoint,
                 coverage = coverage, pop.pred = pop.pred)
}

#' @title Get predictions for new subjects based on a \code{panelTVP.ZINB} object
#'
#' @description
#'  This \code{predict} function simulates data and computes summary statistics based on the
#'   posterior predictive distribution for new data under a Zero-Inflated Negative
#'   Binomial likelihood.
#'
#' @param object an object of class \code{panelTVP.ZINB}
#' @param X_nb.new a matrix or data frame consisting of new data for the same variables
#'  that were used for fitting the model. The first column must contain a 1, when
#'  the fitted model contains an intercept in the Negative Binomial part (count component)
#' @param X_logit.new a matrix or data frame consisting of new data for the same variables
#'  that were used for fitting the model. The first column must contain a 1, when
#'  the fitted model contains an intercept in the Logit part (zero-inflation component)
#' @param timepoint a numeric scalar indicating the time point for which predictions should
#'  be made, i.e., predictions for a given data set are only made for one specific time point.
#'  In case you want predictions for additional time points you need to repeatedly call this
#'  function.
#' @param coverage coverage probability for prediction intervals - defaults to 95 percent coverage
#' @param pop.pred logical value, if TRUE population-based predictions are made,
#'  that ignore the random effects structure, if FALSE the random effects structure
#'  is included as well where the unknown factor scores are sampled from their
#'  standard Normal prior - defaults to FALSE
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @returns a list containing the following elements
#' \itemize{
#'   \item \code{predictive.distribution}: a matrix where each row contains draws
#'    from the posterior predictive distribution for the corresponding observation
#'   \item \code{predictive.summary}: posterior mean as well as HPD interval based on
#'    posterior predictive distribution for each new observation
#' }
#' @exportS3Method predict panelTVP.ZINB
#' @examples
#' # Predictions based on an object of class panelTVP.ZINB
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
#' # setting up design matrix for predicting two new observations
#' X_nb.new <- data.frame(cbind(c(1,1), # column of 1's for the intercept
#'                              c(2,0),
#'                              c(4,3),
#'                              c(0,0)))
#' colnames(X_nb.new) <- colnames(res.zinb$data$X_nb)
#' X_logit.new <- data.frame(cbind(c(1,1), # column of 1's for the intercept
#'                                 c(0.1,-0.5),
#'                                 c(1,2),
#'                                 c(4,3)))
#' colnames(X_logit.new) <- colnames(res.zinb$data$X_logit)
#' # prediction for the 2nd panel wave
#' pp <- predict(res.zinb, X_nb.new = X_nb.new, X_logit.new = X_logit.new, timepoint = 2)
#' plot(density(pp$predictive.distribution[1,]))
#' plot(density(pp$predictive.distribution[2,]))
#' pp$predictive.summary
predict.panelTVP.ZINB <- function(object, X_nb.new, X_logit.new, timepoint,
                                  coverage = 0.95, pop.pred = FALSE, ...){
  check.predict_ZINB(model = object, X_nb.new = X_nb.new, X_logit.new = X_logit.new,
                     timepoint = timepoint, coverage = coverage, pop.pred = pop.pred)
  pred.helper_ZINB(model = object, X_nb.new = X_nb.new, X_logit.new = X_logit.new,
                      timepoint = timepoint, coverage = coverage, pop.pred = pop.pred)
}

pred.helper <- function(model, X.new, timepoint,
                           coverage = 0.95, pop.pred = FALSE){

  if(is.data.frame(X.new)) X.new <- as.matrix(X.new)
  cla <- class(model)
  n <- nrow(X.new)
  d <- ncol(X.new)
  t <- timepoint
  S <- nrow(model$mcmc)
  mcmc.beta <- model$mcmc[, paste0("beta_t", 1:d, t)]
  if(sum(startsWith(colnames(model$mcmc), "lambda_t")) == 1){ # cps
    mcmc.lamb <- model$mcmc[, "lambda_t"]
  } else{
    mcmc.lamb <- model$mcmc[, paste0("lambda_t", t)]
  }
  if(cla == "panelTVP.Gaussian") sigma2 <- model$mcmc[,"sigma2"]
  if(cla == "panelTVP.NegBin") r <- model$mcmc[,"r"]
  y.future <- matrix(nrow = n, ncol = S)
  mcmc.eta <- matrix(nrow = n, ncol = S)
  LO <- numeric(n); UP <- numeric(n)
  if(pop.pred){
    for(s in 1:S){
      mcmc.eta[, s] <- X.new %*% mcmc.beta[s,]
      if(cla == "panelTVP.Gaussian"){
        y.future[, s] <- rnorm(n, mean = mcmc.eta[, s], sd = sqrt(sigma2[s]))
      }
      if(cla == "panelTVP.Probit"){
        y.future[, s] <- rbinom(n, size = 1, prob = pnorm(mcmc.eta[,s]))
      }
      if(cla == "panelTVP.Logit"){
        y.future[, s] <- rbinom(n, size = 1, prob = plogis(mcmc.eta[,s]))
      }
      if(cla == "panelTVP.NegBin"){
        y.future[, s] <- MASS::rnegbin(n, mu = r[s] * exp(mcmc.eta[,s]), theta = r[s])
      }
    }
  } else{
    for(s in 1:S){
      fi <- rnorm(n) # sampling from the prior for new subjects
      mcmc.eta[, s] <- X.new %*% mcmc.beta[s,] + mcmc.lamb[s] * fi
      if(cla == "panelTVP.Gaussian"){
        y.future[, s] <- rnorm(n, mean = mcmc.eta[, s], sd = sqrt(sigma2[s]))
      }
      if(cla == "panelTVP.Probit"){
        y.future[, s] <- rbinom(n, size = 1, prob = pnorm(mcmc.eta[,s]))
      }
      if(cla == "panelTVP.Logit"){
        y.future[, s] <- rbinom(n, size = 1, prob = plogis(mcmc.eta[,s]))
      }
      if(cla == "panelTVP.NegBin"){
        y.future[, s] <- MASS::rnegbin(n, mu = r[s] * exp(mcmc.eta[,s]), theta = r[s])
      }
    }
  }

  # compute summary statistics
  pm <- rowMeans(y.future)
  for(i in 1:n){
    mc <- coda::as.mcmc(y.future[i,])
    hpd <- coda::HPDinterval(mc, prob = coverage)
    LO[i] <- hpd[1]
    UP[i] <- hpd[2]
  }

  res <- list(predictive.distribution = y.future,
              predictive.summary = data.frame(LO = LO, mean = pm, UP = UP))

  return(res)

}

pred.helper_ZINB <- function(model, X_nb.new, X_logit.new, timepoint,
                                coverage = 0.95, pop.pred = FALSE){

  if(is.data.frame(X_nb.new)) X_nb.new <- as.matrix(X_nb.new)
  if(is.data.frame(X_logit.new)) X_logit.new <- as.matrix(X_logit.new)
  n <- nrow(X_nb.new)
  d_nb <- ncol(X_nb.new)
  d_logit <- ncol(X_logit.new)
  t <- timepoint
  S <- nrow(model$mcmc_nb)
  mcmc.beta_nb <- model$mcmc_nb[, paste0("beta_t", 1:d_nb, t)]
  mcmc.beta_logit <- model$mcmc_logit[, paste0("beta_t", 1:d_logit, t)]
  if(sum(startsWith(colnames(model$mcmc_nb), "lambda_t")) == 1){ # cps
    mcmc.lamb_nb <- model$mcmc_nb[, "lambda_t"]
  } else{
    mcmc.lamb_nb <- model$mcmc_nb[, paste0("lambda_t", t)]
  }
  if(sum(startsWith(colnames(model$mcmc_logit), "lambda_t")) == 1){ # cps
    mcmc.lamb_logit <- model$mcmc_logit[, "lambda_t"]
  } else{
    mcmc.lamb_logit <- model$mcmc_logit[, paste0("lambda_t", t)]
  }
  r <- model$mcmc_nb[,"r"]
  y.future <- matrix(nrow = n, ncol = S)
  mcmc.eta_nb <- matrix(nrow = n, ncol = S)
  mcmc.eta_logit <- matrix(nrow = n, ncol = S)
  LO <- numeric(n); UP <- numeric(n)
  if(pop.pred){
    for(s in 1:S){
      mcmc.eta_nb[, s] <- X_nb.new %*% mcmc.beta_nb[s,]
      mcmc.eta_logit[, s] <- X_logit.new %*% mcmc.beta_logit[s,]
      w <- rbinom(n, size = 1, prob = plogis(mcmc.eta_logit[,s]))
      y.future[, s] <- ifelse(w == 1, MASS::rnegbin(n, mu = r[s] * exp(mcmc.eta_nb[,s]),
                                              theta = r[s]), 0)
    }
  } else{
    for(s in 1:S){
      fi_nb <- rnorm(n)
      fi_logit <- rnorm(n)
      mcmc.eta_nb[, s] <- X_nb.new %*% mcmc.beta_nb[s,] +
        mcmc.lamb_nb[s] * fi_nb
      mcmc.eta_logit[, s] <- X_logit.new %*% mcmc.beta_logit[s,] +
        mcmc.lamb_logit[s] * fi_logit
      w <- rbinom(n, size = 1, prob = plogis(mcmc.eta_logit[,s]))
      y.future[, s] <- ifelse(w == 1, MASS::rnegbin(n, mu = r[s] * exp(mcmc.eta_nb[,s]),
                                              theta = r[s]), 0)
    }
  }

  # compute summary statistics
  pm <- rowMeans(y.future)
  for(i in 1:n){
    mc <- coda::as.mcmc(y.future[i,])
    hpd <- coda::HPDinterval(mc, prob = coverage)
    LO[i] <- hpd[1]
    UP[i] <- hpd[2]
  }

  res <- list(predictive.distribution = y.future,
              predictive.summary = data.frame(LO = LO, mean = pm, UP = UP))

  return(res)

}

