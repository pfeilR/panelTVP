#' Fit a Bayesian panel model with time-varying parameters
#'
#' @param formula The usual formula argument in regression methods, e.g., as in [lm()]
#'  (no default)
#' @param data A data.frame containing the variables specified in the formula argument
#'  (no default)
#' @param model A character indicating which model you want to estimate
#'  This parameter is either 'Gaussian', 'Probit', 'Logit' or 'NegBin' depending on
#'  whether you want to fit a model for Gaussian, Probit, Logit or Negative Binomial
#'  response data. In case you want to estimate a Zero-Inflated Negative Binomial
#'  model, please use the function [`panelTVP_ZINB()`], which was designed for exactly
#'  this situation
#' @param prior.reg A list of arguments for estimating the parameters of the regression
#'  part of the model (see details)
#' @param prior.var A list of arguments for estimating the homoscedastic error variance
#'  in a Gaussian/Normal model. For other models, this argument is ignored (see details)
#' @param prior.load A list of arguments for estimating the parameters of the factor
#'  part of the model (see details)
#' @param mcmc.opt A list containing information for the overall sampler (see details)
#' @param settings.NegBin A list containing information for sampling the dispersion
#'  parameter r in the Negative Binomial model. For other response distributions,
#'  this argument is ignored (see details)
#' @param HPD.coverage Coverage probability of highest posterior density intervals
#'  (default yields 95 percent coverage)
#'
#' @returns nix
#' @export
#'
#' @import stats
#'
#' @examples 1
panelTVP <- function(formula,
                     data,
                     model,
                     prior.reg = list(
                       e1 = 0.001, e2 = 0.001, d1 = 0.001, d2 = 0.001,
                       b_tau = 10, nu_tau = 5, b_xi = 10, nu_xi = 5,
                       a_tau = 1, kappa_tau = 10, a_xi = 1, kappa_xi = 10,
                       iota.reg_tau = 1, iota.reg_xi = 1,
                       learn_a_tau = TRUE, learn_a_xi = TRUE,
                       tau.target.rate = 0.44, xi.target.rate = 0.44,
                       learn_kappa_tau = TRUE, learn_kappa_xi = TRUE,
                       type = "rw2", c = 1, B0 = 100000
                     ),
                     prior.var = list(
                       learn.C0.hyp = list(g0 = 5, G0 = 3.333333), c0 = 2.5
                     ),
                     prior.load = list(
                       e1 = 0.001, e2 = 0.001, d1 = 0.001, d2 = 0.001,
                       b_phi = 10, nu_phi = 5, b_zeta = 10, nu_zeta = 5,
                       a_phi = 1, kappa_phi = 10, a_zeta = 1, kappa_zeta = 10,
                       learn_kappa_phi = TRUE, learn_kappa_zeta = TRUE,
                       type = "rw2", c = 1, L0 = 1
                     ),
                     mcmc.opt = list(
                       chain.length = 12000, burnin = 2000, thin = 10, asis = TRUE
                     ),
                     settings.NegBin = list(
                       iota.r = 0.1581139, r.a = 0.1, r.b = 10, r.target.rate = 0.44,
                       slice = TRUE, r.alpha = 2, r.beta = 1, expansion.steps = 20,
                       width = 1, p.overrelax = 0, accuracy.overrelax = 20
                     ),
                     HPD.coverage = 0.95
){

  # Initialization -------------------------------------------------------------

  if(prior.load$type == "cps"){
    tv.load = FALSE
  } else{
    tv.load = TRUE
  }

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
                namestau, namesxi, c("a_tau","kappa_tau","a_xi","kappa_xi"),
                namessgma2, nameslambdat)
  }
  if(prior.load$type=="rw1" | prior.load$type=="rw2"){
    cnames <- c(cnames,"lambda","psi","phi2","zeta2", "a_phi", "kappa_phi", "a_zeta", "kappa_zeta")
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
  hyperpara <- c("a_xi", "a_tau", "kappa_xi", "kappa_tau", "kappa_zeta", "kappa_phi")
  part <- c(rep("regression part", 4), rep("factor part", 2))
  learn <- c(prior.reg$learn_a_xi, prior.reg$learn_a_tau,
             prior.reg$learn_kappa_xi, prior.reg$learn_kappa_tau,
             prior.load$learn_kappa_zeta, prior.load$learn_kappa_phi)
  if(!(prior.reg$type %in% c("rw1", "rw2"))) learn[1:4] <- NA
  if(!(prior.load$type %in% c("rw1", "rw2"))) learn[5:6] <- NA
  result$learning.settings <- cbind(hyperpara, part, learn)
  colnames(result$learning.settings) <- c("hyperparameter", "model.part", "learned?")

  # adding mcmc setting to output (incl. ASIS Boolean)
  result$mcmc.settings <- mcmc.opt

  return(result)

}
