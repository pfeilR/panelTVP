#' @title Summary output for a \code{panelTVP.Gaussian} object
#'
#' @description
#'   This \code{summary} function prints out a table that contains HPD-intervals,
#'   posterior means and standard deviations based on the posterior distributions of the
#'   time-varying parameters, i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T}. The results are sorted either by covariate
#'   (default) or by time point.
#'
#' @param object an object of class \code{panelTVP.Gaussian}
#' @param by a single character that is either 'timepoint' or 'covariate'
#'  to sort the output either by time point or covariate, respectively (default is covariate)
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method summary panelTVP.Gaussian
#' @examples
#' # Computing summary table for object of class panelTVP.Gaussian
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
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#' summary(res.gaussian) # default = ordering by covariate
#' summary(res.gaussian, by = "timepoint") # ordering by time point
summary.panelTVP.Gaussian <- function(object, by = "covariate", ...){
  if(length(by)>1 | !is.character(by)){
    stop("by is a scalar character")
  }
  if((!by %in% c("timepoint", "covariate"))){
    stop("you can only sort it by timepoint or covariate")
  }
  cat("\n------------------------------------------------------------------------------
Posterior Summary of the Bayesian Normal Model with Time-Varying Coefficients:
------------------------------------------------------------------------------\n")
  cat(craft.summary(object, by = by))
}

#' @title Summary output for a \code{panelTVP.Probit} object
#'
#' @description
#'   This \code{summary} function prints out a table that contains HPD-intervals,
#'   posterior means and standard deviations based on the posterior distributions of the
#'   time-varying parameters, i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T}. The results are sorted either by covariate
#'   (default) or by time point.
#'
#' @param object an object of class \code{panelTVP.Probit}
#' @param by a single character that is either 'timepoint' or 'covariate'
#'  to sort the output either by time point or covariate, respectively (default is covariate)
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method summary panelTVP.Probit
#' @examples
#' # Computing summary table for object of class panelTVP.Probit
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
#'                        id = sim.probit$observed$id,
#'                        t = sim.probit$observed$t,
#'                        mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                        model = "Probit")
#' summary(res.probit) # default = ordering by covariate
#' summary(res.probit, by = "timepoint") # ordering by time point
summary.panelTVP.Probit <- function(object, by = "covariate", ...){
  if(length(by)>1 | !is.character(by)){
    stop("by is a scalar character")
  }
  if((!by %in% c("timepoint", "covariate"))){
    stop("you can only sort it by timepoint or covariate")
  }
  cat("------------------------------------------------------------------------------
Posterior Summary of the Bayesian Probit Model with Time-Varying Coefficients:
------------------------------------------------------------------------------\n")
  cat(craft.summary(object, by = by))
}

#' @title Summary output for a \code{panelTVP.Logit} object
#'
#' @description
#'   This \code{summary} function prints out a table that contains HPD-intervals,
#'   posterior means and standard deviations based on the posterior distributions of the
#'   time-varying parameters, i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T}. The results are sorted either by covariate
#'   (default) or by time point.
#'
#' @param object an object of class \code{panelTVP.Logit}
#' @param by a single character that is either 'timepoint' or 'covariate'
#'  to sort the output either by time point or covariate, respectively (default is covariate)
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method summary panelTVP.Logit
#' @examples
#' # Computing summary table for object of class panelTVP.Logit
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
#'                       id = sim.logit$observed$id,
#'                       t = sim.logit$observed$t,
#'                       mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                       model = "Logit")
#' summary(res.logit) # default = ordering by covariate
#' summary(res.logit, by = "timepoint") # ordering by time point
summary.panelTVP.Logit <- function(object, by = "covariate", ...){
  if(length(by)>1 | !is.character(by)){
    stop("by is a scalar character")
  }
  if((!by %in% c("timepoint", "covariate"))){
    stop("you can only sort it by timepoint or covariate")
  }
  cat("------------------------------------------------------------------------------
Posterior Summary of the Bayesian Logit Model with Time-Varying Coefficients:
-----------------------------------------------------------------------------\n")
  cat(craft.summary(object, by = by))
}

#' @title Summary output for a \code{panelTVP.NegBin} object
#'
#' @description
#'   This \code{summary} function prints out a table that contains HPD-intervals,
#'   posterior means and standard deviations based on the posterior distributions of the
#'   time-varying parameters, i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T}. The results are sorted either by covariate
#'   (default) or by time point.
#'
#' @param object an object of class \code{panelTVP.NegBin}
#' @param by a single character that is either 'timepoint' or 'covariate'
#'  to sort the output either by time point or covariate, respectively (default is covariate)
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method summary panelTVP.NegBin
#' @examples
#' # Computing summary table for object of class panelTVP.NegBin
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
#'                        id = sim.negbin$observed$id,
#'                        t = sim.negbin$observed$t,
#'                        mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                        model = "NegBin")
#' summary(res.negbin) # default = ordering by covariate
#' summary(res.negbin, by = "timepoint") # ordering by time point
summary.panelTVP.NegBin <- function(object, by = "covariate", ...){
  if(length(by)>1 | !is.character(by)){
    stop("by is a scalar character")
  }
  if((!by %in% c("timepoint", "covariate"))){
    stop("you can only sort it by timepoint or covariate")
  }
  cat("-----------------------------------------------------------------------------------------
Posterior Summary of the Bayesian Negative Binomial Model with Time-Varying Coefficients:
-----------------------------------------------------------------------------------------\n")
  cat(craft.summary(object, by = by))
}

#' @title Summary output for a \code{panelTVP.ZINB} object
#'
#' @description
#'   This \code{summary} function prints out a table that contains HPD-intervals,
#'   posterior means and standard deviations based on the posterior distributions of the
#'   time-varying parameters, i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T}. The results are sorted either by covariate (default)
#'   or by time point.
#'
#' @param object an object of class \code{panelTVP.ZINB}
#' @param by a single character that is either 'timepoint' or 'covariate'
#'  to sort the output either by time point or covariate, respectively (default is covariate)
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method summary panelTVP.ZINB
#' @examples
#' # Computing summary table for object of class panelTVP.ZINB
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
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
#' summary(res.zinb) # default = ordering by covariate
#' summary(res.zinb, by = "timepoint") # ordering by time point
summary.panelTVP.ZINB <- function(object, by = "covariate", ...){
  if(length(by)>1 | !is.character(by)){
    stop("by is a scalar character")
  }
  if((!by %in% c("timepoint", "covariate"))){
    stop("you can only sort it by timepoint or covariate")
  }
  cat("-------------------------------------------------------------------------------------------------------
Posterior Summary of the Bayesian Zero-Inflated Negative Binomial Model with Time-Varying Coefficients:
-------------------------------------------------------------------------------------------------------\n")
  if(by == "timepoint"){
    cat(craft.summary_zinb(object))
  } else{
    cat(craft.summary_zinb_by_covariate(object))
  }

}

craft.summary <- function(x, by = by){

  res <- crafti(X = x$data$X, posterior = x$posterior, by = by, ntime = x$data$Tmax)
  nami <- names(res)
  output_lines <- c()
  for (i in 1:length(res)) {
    mod <- round(res[[i]], 4)

    # Header
    if(by == "timepoint"){
      output_lines <- c(
        output_lines,
        strrep("=", 50),
        center_text(paste("Estimates for Timepoint", i)),
        strrep("=", 50),
        ""
      )
    } else{
      output_lines <- c(
        output_lines,
        strrep("=", 50),
        center_text(paste("Estimates for", nami[i])),
        strrep("=", 50),
        ""
      )
    }

    # Model output
    output_lines <- c(output_lines, utils::capture.output(print(mod)), "")

  }

  # Combine into a single string if needed
  output_string <- paste(output_lines, collapse = "\n")
  return(output_string)
}

craft.summary_zinb <- function(x){

  res_nb <- crafti(X = x$data$X_nb, posterior = x$posterior_nb, by = "timepoint",
                   ntime = x$data$Tmax)
  res_logit <- crafti(X = x$data$X_logit, posterior = x$posterior_logit, by = "timepoint",
                      ntime = x$data$Tmax)

  output_lines <- c()
  for (i in 1:length(res_nb)) {
    nb <- round(res_nb[[i]], 4)
    logit <- round(res_logit[[i]], 4)

    # Header
    output_lines <- c(
      output_lines,
      strrep("=", 50),
      center_text(paste("Estimates: Timepoint", i)),
      strrep("=", 50),
      ""
    )

    # Negative Binomial model
    output_lines <- c(output_lines, center_text("---- Negative Binomial Model ----"), "")
    output_lines <- c(output_lines, utils::capture.output(print(nb)), "")

    # Logit model
    output_lines <- c(output_lines, center_text("---- Logit Model ----"), "")
    output_lines <- c(output_lines, utils::capture.output(print(logit)), "", "")
  }

  # Combine into a single string if needed
  output_string <- paste(output_lines, collapse = "\n")
  return(output_string)

}

craft.summary_zinb_by_covariate <- function(x){

  res_nb <- crafti(X = x$data$X_nb, posterior = x$posterior_nb, by = "covariate",
                   ntime = x$data$Tmax)
  res_logit <- crafti(X = x$data$X_logit, posterior = x$posterior_logit, by = "covariate",
                      ntime = x$data$Tmax)
  nami_nb <- names(res_nb)
  nami_logit <- names(res_logit)

  output_lines <- c()
  for (i in 1:length(res_nb)) {
    nb <- round(res_nb[[i]], 4)

    # Header
    output_lines <- c(
      output_lines,
      strrep("=", 50),
      center_text(paste("Estimates for", nami_nb[i])),
      strrep("=", 50),
      ""
    )
    output_lines <- c(output_lines, utils::capture.output(print(nb)), "")

  }
  output_lines <- c(strrep("-", 50), center_text("Negative Binomial Model"), strrep("-", 50),
                    output_lines, strrep("-", 50), center_text("Logit Model"), strrep("-", 50))

  for (i in 1:length(res_logit)) {
    logit <- round(res_logit[[i]], 4)

    # Header
    output_lines <- c(
      output_lines,
      strrep("=", 50),
      center_text(paste("Estimates for", nami_logit[i])),
      strrep("=", 50),
      ""
    )
    output_lines <- c(output_lines, utils::capture.output(print(logit)), "")

  }

  # Combine into a single string if needed
  output_string <- paste(output_lines, collapse = "\n")
  return(output_string)

}

crafti <- function(X, posterior, by = NULL, ntime){

  d <- posterior[,c(1,2,4,5)]
  colnames(d) <- c("Lower (HPD)", "Posterior Mean", "Upper (HPD)", "SD")
  d <- d[startsWith(rownames(d), "beta_t")|startsWith(rownames(d), "lambda_t"),]
  ncov <- sum(startsWith(rownames(d), "beta_t"))/ntime
  d <- as.data.frame(d)
  if(nrow(d[startsWith(rownames(d), "lambda_t"),]) == 1){ # cps prior
    mat <- as.matrix(d[startsWith(rownames(d), "lambda_t"),])
    mat <- mat[rep(1,ntime),]
    rownames(mat) <- paste0("lambda_t", 1:ntime)
    d <- rbind(d[!startsWith(rownames(d), "lambda_t"),], mat)
  }
  d$time <- rep(1:ntime, times = ncov + 1)
  d$cov <- c(rep(1:ncov, each = ntime), rep(ncov+1, ntime))
  if(by == "timepoint"){
    d <- d[order(d$time,d$cov),]
  } else{
    d <- d[order(d$cov,d$time),]
  }
  namesbeta <- colnames(X)
  if(by == "timepoint"){
    cnames <- rep(c(namesbeta, "\u03BB"), times = ntime)
  } else{
    cnames <- rep(paste0("t=",1:ntime), times = ncov+1) # +1 for factor loading
  }
  d <- as.matrix(d)
  rownames(d) <- cnames
  res <- list()
  if(by == "timepoint"){
    for(t in 1:ntime){
      res[[t]] <- d[d[,"time"] == t, 1:4]
      names(res)[t] <- paste("Regression Effects and Factor Loading at Time", t)
    }
  } else{
    covnames <- c(namesbeta, "\u03BB")
    for(i in 1:(ncov+1)){ # +1 for factor loading
      res[[i]] <- d[d[,"cov"] == i, 1:4]
      names(res)[i] <- covnames[i]
    }
  }
  return(res)

}

center_text <- function(text, width = 50) {
  if(nchar(text) > width) {
    warning("At least one variable name is too long for fully displaying it.")
    text <- substr(text, 1, width)
  }
  padding <- floor((width - nchar(text)) / 2)
  paste0(strrep(" ", padding), text)
}
