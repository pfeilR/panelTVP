#' Fit a Bayesian panel model with time-varying parameters for overdispersed and zero-inflated count data
#'
#' @param formula A formula argument consisting of two parts: the count model and
#'  the zero-inflation model (no default, see details)
#' @param data A data.frame containing the variables specified in the formula argument
#'  (no default)
#' @param prior.reg_nb A list of arguments for estimating the parameters of the regression
#'  part of the count model (see details)
#' @param prior.load_nb A list of arguments for estimating the parameters of the factor
#'  part of the count model (see details)
#' @param prior.reg_logit A list of arguments for estimating the parameters of the regression
#'  part of the zero-inflation model (see details)
#' @param prior.load_logit A list of arguments for estimating the parameters of the factor
#'  part of the zero-inflation model (see details)
#' @param mcmc.opt A list containing information for the overall sampler (see details)
#' @param settings.NegBin A list containing information for sampling the dispersion
#'  parameter r of the count model (see details)
#' @param HPD.coverage Coverage probability of highest posterior density intervals
#'  (default yields 95 percent coverage)
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
#'    \item{mcmc.settings}{details on general MCMC sampler}
#'  }
#' @export
#'
#' @examples 0
panelTVP_ZINB <- function(formula,
                          data,
                          prior.reg_nb = list(
                            e1 = 0.001, e2 = 0.001, d1 = 0.001, d2 = 0.001,
                            b_tau = 10, nu_tau = 5, b_xi = 10, nu_xi = 5,
                            a_tau = 1, kappa_tau = 10, a_xi = 1, kappa_xi = 10,
                            iota.reg_tau = 1, iota.reg_xi = 1,
                            learn_a_tau = TRUE, learn_a_xi = TRUE,
                            tau.target.rate = 0.44, xi.target.rate = 0.44,
                            learn_kappa_tau = TRUE, learn_kappa_xi = TRUE,
                            type = "rw2", c = 1, B0 = 100000
                          ),
                          prior.load_nb = list(
                            e1 = 0.001, e2 = 0.001, d1 = 0.001, d2 = 0.001,
                            b_phi = 10, nu_phi = 5, b_zeta = 10, nu_zeta = 5,
                            a_phi = 1, kappa_phi = 10, a_zeta = 1, kappa_zeta = 10,
                            learn_kappa_phi = TRUE, learn_kappa_zeta = TRUE,
                            type = "rw2", c = 1, L0 = 1
                          ),
                          prior.reg_logit = list(
                            e1 = 0.001, e2 = 0.001, d1 = 0.001, d2 = 0.001,
                            b_tau = 10, nu_tau = 5, b_xi = 10, nu_xi = 5,
                            a_tau = 1, kappa_tau = 10, a_xi = 1, kappa_xi = 10,
                            iota.reg_tau = 1, iota.reg_xi = 1,
                            learn_a_tau = TRUE, learn_a_xi = TRUE,
                            tau.target.rate = 0.44, xi.target.rate = 0.44,
                            learn_kappa_tau = TRUE, learn_kappa_xi = TRUE,
                            type = "rw2", c = 1, B0 = 100000
                          ),
                          prior.load_logit = list(
                            e1 = 0.001, e2 = 0.001, d1 = 0.001, d2 = 0.001,
                            b_phi = 10, nu_phi = 5, b_zeta = 10, nu_zeta = 5,
                            a_phi = 1, kappa_phi = 10, a_zeta = 1, kappa_zeta = 10,
                            learn_kappa_phi = TRUE, learn_kappa_zeta = TRUE,
                            type = "rw2", c = 1, L0 = 1
                          ),
                          mcmc.opt = list(
                            chain.length = 12000, thin = 10, burnin = 2000, asis = TRUE
                          ),
                          settings.NegBin = list(
                            iota.r = 0.1581139, r.a = 0.1, r.b = 10, r.target.rate = 0.44,
                            slice = TRUE, r.alpha = 2, r.beta = 1, expansion.steps = 20,
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
                   namestau_nb, namesxi_nb, c("a_tau","kappa_tau","a_xi","kappa_xi"),
                   nameslambdat_nb)
  }
  if(prior.load_nb$type=="rw1" | prior.load_nb$type=="rw2"){
    cnames_nb <- c(cnames_nb,"lambda","psi","phi2","zeta2", "a_phi", "kappa_phi",
                   "a_zeta", "kappa_zeta")
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
                      namestau_logit, namesxi_logit, c("a_tau","kappa_tau","a_xi","kappa_xi"),
                      nameslambdat_logit)
  }
  col_res_logit <- length(cnames_logit)
  if(prior.load_logit$type=="rw1" | prior.load_logit$type=="rw2"){
    cnames_logit <- c(cnames_logit,"lambda","psi","phi2","zeta2", "a_phi", "kappa_phi",
                      "a_zeta", "kappa_zeta")
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
  hyperpara <- c("a_xi", "a_tau", "kappa_xi", "kappa_tau", "kappa_zeta", "kappa_phi")
  part <- c(rep("regression part", 4), rep("factor part", 2))
  learn_nb <- c(prior.reg_nb$learn_a_xi, prior.reg_nb$learn_a_tau,
                prior.reg_nb$learn_kappa_xi, prior.reg_nb$learn_kappa_tau,
                prior.load_nb$learn_kappa_zeta, prior.load_nb$learn_kappa_phi)
  if(!(prior.reg_nb$type %in% c("rw1", "rw2"))) learn_nb[1:4] <- NA
  if(!(prior.load_nb$type %in% c("rw1", "rw2"))) learn_nb[5:6] <- NA
  learn_logit <- c(prior.reg_logit$learn_a_xi, prior.reg_logit$learn_a_tau,
                   prior.reg_logit$learn_kappa_xi, prior.reg_logit$learn_kappa_tau,
                   prior.load_logit$learn_kappa_zeta, prior.load_logit$learn_kappa_phi)
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

