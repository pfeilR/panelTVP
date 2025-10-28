# check input arguments of sim_panelTVP ----------------------------------------

check_sim <- function(n,
                      Tmax,
                      model,
                      beta = NULL,
                      theta = NULL,
                      lambda = NULL,
                      psi = NULL,
                      r = NULL,
                      sigma2 = NULL,
                      beta.nb = NULL,
                      theta.nb = NULL,
                      lambda.nb = NULL,
                      psi.nb = NULL,
                      beta.logit = NULL,
                      theta.logit = NULL,
                      lambda.logit = NULL,
                      psi.logit = NULL){

  if(is.null(n) || !is.numeric(n) || length(n) != 1 || n %% 1 != 0 || !is.finite(n)){
    stop("Argument 'n' must be a single value that represents the number of subjects")
  }
  if(is.null(Tmax) || !is.numeric(Tmax) || length(Tmax) != 1 || Tmax %% 1 != 0 || !is.finite(Tmax) || Tmax <= 2){
    stop("Argument 'Tmax' must be a single value > 2 that represents the number of repeated measurements")
  }
  if(is.null(model) || !(model %in% c("Gaussian", "Probit", "Logit", "NegBin", "ZINB"))){
    stop("Argument 'model' must be either 'Gaussian', 'Probit', 'Logit', 'NegBin' or 'ZINB'")
  }
  if(model %in% c("NegBin", "ZINB") && is.null(r)){
    stop("Argument 'r' must be specified when model is either 'NegBin' or 'ZINB'")
  }
  if(model == "Gaussian" && is.null(sigma2)){
    stop("Argument 'sigma2' must be specified when model is 'Gaussian'")
  }
  if(model %in% c("Gaussian", "Probit", "Logit", "NegBin") &&
     any(sapply(list(beta, theta, lambda, psi), is.null))){
    stop("When model is either 'Gaussian', 'Probit', 'Logit', or 'NegBin' arguments 'beta', 'theta', 'lambda' and 'psi' have to be specified")
  }
  if(model == "ZINB" && any(sapply(list(beta.nb, theta.nb, lambda.nb, psi.nb,
                                        beta.logit, theta.logit, lambda.logit, psi.logit), is.null))){
    stop("When model is 'ZINB' arguments 'beta.nb', 'theta.nb', 'lambda.nb', 'psi.nb', 'beta.logit', 'theta.logit', 'lambda.logit', 'psi.logit' have to be specified")
  }
  if(all(sapply(list(beta, theta, lambda, psi), function(x) !is.null(x)))){
    if(any(sapply(list(beta, theta, lambda, psi), function(x) length(x) == 0))){
      stop("Arguments 'beta', 'theta', 'lambda', 'psi' must be of positive length")
    }
    if(length(beta) != length(theta)) stop("Argument 'beta' must be of same length as argument 'theta'")
    if(!is.numeric(beta) || any(!is.finite(beta)) || length(beta) < 2) stop("Argument 'beta' must be vector of finite values with 2 or more elements")
    if(!is.numeric(theta) || any(!is.finite(theta)) || length(theta) < 2) stop("Argument 'theta' must be vector of finite values with 2 or more elements")
    if(!is.numeric(lambda) || length(lambda) > 1 || !is.finite(lambda)) stop("Argument 'lambda' must be finite numeric")
    if(!is.numeric(psi) || length(psi) > 1 || !is.finite(psi)) stop("Argument 'psi' must be finite numeric")
  }
  if(all(sapply(list(beta.nb, theta.nb, lambda.nb, psi.nb,
                     beta.logit, theta.logit, lambda.logit, psi.logit), function(x) !is.null(x)))){
    if(any(sapply(list(beta.nb, theta.nb, lambda.nb, psi.nb,
                       beta.logit, theta.logit, lambda.logit, psi.logit), function(x) length(x) == 0))){
      stop("Arguments 'beta.nb', 'theta.nb', 'lambda.nb', 'psi.nb', 'beta.logit', 'theta.logit', 'lambda.logit', 'psi.logit' must be of positive length")
    }
    if(length(beta.nb) != length(theta.nb)) stop("Argument 'beta.nb' must be of same length as argument 'theta.nb'")
    if(!is.numeric(beta.nb) || any(!is.finite(beta.nb)) || length(beta.nb) < 2) stop("Argument 'beta.nb' must be vector of finite values with 2 or more elements")
    if(!is.numeric(theta.nb) || any(!is.finite(theta.nb)) ||length(theta.nb) < 2) stop("Argument 'theta.nb' must be vector of finite values with 2 or more elements")
    if(!is.numeric(lambda.nb) || length(lambda.nb) > 1 || !is.finite(lambda.nb)) stop("Argument 'lambda.nb' must be finite numeric")
    if(!is.numeric(psi.nb) || length(psi.nb) > 1 || !is.finite(psi.nb)) stop("Argument 'psi.nb' must be finite numeric")
    if(length(beta.logit) != length(theta.logit)) stop("Argument 'beta.logit' must be of same length as argument 'theta.logit'")
    if(!is.numeric(beta.logit) || any(!is.finite(beta.logit)) ||length(beta.logit) < 2) stop("Argument 'beta.logit' must be vector of finite values with 2 or more elements")
    if(!is.numeric(theta.logit) || any(!is.finite(theta.logit)) ||length(theta.logit) < 2) stop("Argument 'theta.logit' must be vector of finite values with 2 or more elements")
    if(!is.numeric(lambda.logit) || length(lambda.logit) > 1 || !is.finite(lambda.logit)) stop("Argument 'lambda.logit' must be finite numeric")
    if(!is.numeric(psi.logit) || length(psi.logit) > 1 || !is.finite(psi.logit)) stop("Argument 'psi.logit' must be finite numeric")
  }
  if(!is.null(r) && (!is.numeric(r) || length(r) > 1 || r < 0 || !is.finite(r))){
    stop("Argument 'r' must be a finite, positive numeric")
  }
  if(!is.null(sigma2) && (!is.numeric(sigma2) || length(sigma2) > 1 || sigma2 < 0 || !is.finite(sigma2))){
    stop("Argument 'sigma2' must be a finite, positive numeric")
  }

}

# check input arguments of panelTVP --------------------------------------------

check.panelTVP <- function(formula, data, id, t, model, prior.reg, prior.var, prior.load,
                           prior.reg_nb, prior.load_nb, prior.reg_logit, prior.load_logit,
                           mcmc.opt, settings.NegBin, HPD.coverage, R.WAIC, posterior.predictive.matrix,
                           random.effects, progress.bar){

  if(is.null(formula) || is.null(data)) stop("Arguments 'formula' and 'data' are required with no defaults.")

  # checking model type
  if(is.null(model) || !model %in% c("Gaussian", "Probit", "Logit", "NegBin", "ZINB")){
    stop("Argument 'model' must be either 'Gaussian', 'Probit', 'Logit', 'NegBin' or 'ZINB'.")
  }
  # checking formula (together with model)
  if(!inherits(formula, "formula")) {
    stop("Argument 'formula' must be a formula object.")
  }
  if(model != "ZINB" && any(grepl("\\|", deparse(formula)))){
    stop("The character '|' in argument 'formula' is only valid for Zero-Inflated Negative Binomial regression.")
  }
  if(model == "ZINB" && sum(grepl("\\|", deparse(formula))) != 1){
    stop("Argument 'formula' needs exactly one '|' character for Zero-Inflated Negative Binomial regression.")
  }
  # data, id, t
  if(is.null(data) || !is.data.frame(data)){
    stop("Argument 'data' must be a (Tn x d) data frame.")
  }
  if(model != "ZINB"){
    vars <- all.vars(formula)
    if(sum(!(vars %in% colnames(data))) > 0){
      if(formula[[3]] != "."){
        stop("There are variables in your 'formula' argument that are not contained in 'data'.")
      }
    }
  } else{
    counti <- formula[[3]][[2]]
    vars.counti <- all.vars(counti)
    if(sum(!(vars.counti %in% colnames(data))) > 0){
      if(counti != "."){
        stop("There are variables in the count part of your 'formula' argument that are not contained in 'data'.")
      }
    }
    logi <- formula[[3]][[3]]
    vars.logi <- all.vars(logi)
    if(sum(!(vars.logi %in% colnames(data))) > 0){
      if(logi != "."){
        stop("There are variables in the zero-inflation part of your 'formula' argument that are not contained in 'data'.")
      }
    }
  }
  vars <- all.vars(formula)
  if("t" %in% vars) stop("Covariates in the formula are not allowed to be called 't'. This is for the time index. Please rename the variable.")
  if("id" %in% vars) stop("Covariates in the formula are not allowed to be called 'id'. This is for the subject index. Please rename the variable.")
  if(is.null(t) && (length(t) != nrow(data) || !is.numeric(t) || sum(t %% 1) != 0 || sum(!is.finite(t)) != 0)){
    stop("Argument 't' must be an integer-valued vector with length equal to the number of observations in 'data'.")
  }
  if(is.null(id) && (length(id) != nrow(data) || !is.numeric(id) || sum(id %% 1) != 0 || sum(!is.finite(id)) != 0)){
    stop("Argument 'id' must be an integer-valued vector with length equal to the number of observations in 'data'.")
  }

  # check if panel is balanced
  if(length(t) %% max(t) != 0){
    stop("Panel must be balanced, i.e., the same number of repeated measurements for every subject.")
  }

  # response variable check
  resp <- data[, as.character(formula[[2]])]
  if(length(unique(resp)) == 1){
    stop("There is no variation in your response variable.")
  }
  if(model %in% c("Probit", "Logit")){
    if(!is.numeric(resp) || !all(resp %in% c(0,1))){
      stop("For binary regression, response must be a numeric with categories coded as 0 and 1.")
    }
  }
  if(model %in% c("NegBin", "ZINB")){
    if(!is.numeric(resp)){
      stop("Response must be numeric for count data regression.")
    }
    if(any(resp < 0)){
      stop("Response must be positive for count data regression.")
    }
    if(any(abs(resp - round(resp)) > .Machine$double.eps^0.5)){
      stop("Response must be integer-valued (count data) for count data regression.")
    }
  }

  # hyperparameter checks for non-ZINB models
  if(model != "ZINB"){

    if(as.character(formula[[2]]) %in% as.character(formula[[3]])){
      stop("Response variable is also contained as predictor.")
    }

    check.prior.reg(prior.reg)
    if(model == "Gaussian"){
      # prior.var.check
      check.hyper_double.positive(prior.var$learn.C0.hyp$g0, "g0 in prior.var")
      check.hyper_double.positive(prior.var$learn.C0.hyp$G0, "G0 in prior.var")
      check.hyper_double.positive(prior.var$c0, "c0 in prior.var")
    }
    check.prior.load(prior.load)

  } else{ # ZINB-Checks

    if(as.character(formula[[2]]) %in% c(as.character(formula[[3]][[2]]),
                                         as.character(formula[[3]][[3]]))){
      stop("Response variable is also contained as predictor.")
    }

    check.prior.reg(prior.reg_nb)
    check.prior.load(prior.load_nb)
    check.prior.reg(prior.reg_logit)
    check.prior.load(prior.load_logit)

  }

  # MCMC setting checks
  check.mcmc.opt(mcmc.opt)

  # settings.NegBin
  if(model %in% c("NegBin", "ZINB")) check.settings.NegBin(settings.NegBin)

  # HPD coverage
  if(is.null(HPD.coverage) ||!is.numeric(HPD.coverage) || HPD.coverage < 0 || HPD.coverage > 1)
    stop("Argument HPD.coverage must be valid probability.")

  # R.WAIC
  if(is.null(R.WAIC) ||!is.numeric(R.WAIC) || length(R.WAIC) != 1 || R.WAIC %% 1 != 0 ||
     !is.finite(R.WAIC) || R.WAIC < 2)
    stop("Argument 'R.WAIC' must be a single, positive integer and has to be at least 2.")

  # posterior.predictive.matrix
  if(is.null(posterior.predictive.matrix) || !is.logical(posterior.predictive.matrix) ||
     length(posterior.predictive.matrix) != 1 || !is.finite(posterior.predictive.matrix)){
    stop("Argument 'posterior.predictive.matrix' must be a single logical value.")
  }

  # random.effects
  if(is.null(random.effects) || !is.logical(random.effects) ||
     length(random.effects) != 1 || !is.finite(random.effects)){
    stop("Argument 'random.effects' must be a single logical value.")
  }

  # progress.bar
  if(is.null(progress.bar) || !is.logical(progress.bar) ||
     length(progress.bar) != 1 || !is.finite(progress.bar)){
    stop("Argument 'progress.bar' must be a single logical value.")
  }


}

# helper for check.panelTVP ----------------------------------------------------
check.hyper_double.positive <- function(x, argument){
  if(is.null(x) || !is.numeric(x) || length(x) != 1 || !is.finite(x) || x <= 0)
    stop(paste("Argument", argument, "must be a positive and finite numeric scalar."))
}

check.prior.reg <- function(prior.reg){
  if(is.null(prior.reg)) stop("Your regression prior setting is not allowed to be NULL.")
  check.hyper_double.positive(prior.reg$a.tau, "a.tau in regression prior")
  check.hyper_double.positive(prior.reg$a.xi, "a.xi in regression prior")
  check.hyper_double.positive(prior.reg$alpha.a.tau, "alpha.a.tau in regression prior")
  check.hyper_double.positive(prior.reg$alpha.a.xi, "alpha.a.xi in regression prior")
  check.hyper_double.positive(prior.reg$beta.a.tau, "beta.a.tau in regression prior")
  check.hyper_double.positive(prior.reg$beta.a.xi, "beta.a.xi in regression prior")
  check.hyper_double.positive(prior.reg$iota.a.tau, "iota.a.tau in regression prior")
  check.hyper_double.positive(prior.reg$iota.a.xi, "iota.a.xi in regression prior")
  check.hyper_double.positive(prior.reg$c.tau, "c.tau in regression prior")
  check.hyper_double.positive(prior.reg$c.xi, "c.xi in regression prior")
  check.hyper_double.positive(prior.reg$alpha.c.tau, "alpha.c.tau in regression prior")
  check.hyper_double.positive(prior.reg$alpha.c.xi, "alpha.c.xi in regression prior")
  check.hyper_double.positive(prior.reg$beta.c.tau, "beta.c.tau in regression prior")
  check.hyper_double.positive(prior.reg$beta.c.xi, "beta.c.xi in regression prior")
  check.hyper_double.positive(prior.reg$iota.c.tau, "iota.c.tau in regression prior")
  check.hyper_double.positive(prior.reg$iota.c.xi, "iota.c.xi in regression prior")
  check.hyper_double.positive(prior.reg$kappa.tau, "kappa.tau in regression prior")
  check.hyper_double.positive(prior.reg$kappa.xi, "kappa.xi in regression prior")
  check.hyper_double.positive(prior.reg$d.tau, "d.tau in regression prior")
  check.hyper_double.positive(prior.reg$d.xi, "d.xi in regression prior")
  check.hyper_double.positive(prior.reg$e.tau, "e.tau in regression prior")
  check.hyper_double.positive(prior.reg$e.xi, "e.xi in regression prior")
  if(!is.logical(prior.reg$learn.a.tau) || length(prior.reg$learn.a.tau) > 1)
    stop("Argument learn.a.tau in regression prior must be a logical scalar.")
  if(!is.logical(prior.reg$learn.a.xi) || length(prior.reg$learn.a.xi) > 1)
    stop("Argument learn.a.xi in regression prior must be a logical scalar.")
  if(!is.logical(prior.reg$learn.c.tau) || length(prior.reg$learn.c.tau) > 1)
    stop("Argument learn.c.tau in regression prior must be a logical scalar.")
  if(!is.logical(prior.reg$learn.c.xi) || length(prior.reg$learn.c.xi) > 1)
    stop("Argument learn.c.xi in regression prior must be a logical scalar.")
  if(prior.reg$learn.a.tau){
    if(is.null(prior.reg$target.rate.a.tau) ||!is.numeric(prior.reg$target.rate.a.tau) ||
       prior.reg$target.rate.a.tau < 0 || prior.reg$target.rate.a.tau > 1)
      stop("Argument target.rate.a.tau in regression prior must be valid probability when learning a.tau.")
  }
  if(prior.reg$learn.a.xi){
    if(is.null(prior.reg$target.rate.a.xi) || !is.numeric(prior.reg$target.rate.a.xi) ||
       prior.reg$target.rate.a.xi < 0 || prior.reg$target.rate.a.xi > 1)
      stop("Argument target.rate.a.xi in regression prior must be valid probability when learning a.xi.")
  }
  if(prior.reg$learn.c.tau){
    if(is.null(prior.reg$target.rate.c.tau) ||!is.numeric(prior.reg$target.rate.c.tau) ||
       prior.reg$target.rate.c.tau < 0 || prior.reg$target.rate.c.tau > 1)
      stop("Argument target.rate.c.tau in regression prior must be valid probability when learning c.tau.")
  }
  if(prior.reg$learn.c.xi){
    if(is.null(prior.reg$target.rate.c.xi) || !is.numeric(prior.reg$target.rate.c.xi) ||
       prior.reg$target.rate.c.xi < 0 || prior.reg$target.rate.c.xi > 1)
      stop("Argument target.rate.c.xi in regression prior must be valid probability when learning c.xi.")
  }
  if(is.null(prior.reg$learn.kappa.tau) || !is.logical(prior.reg$learn.kappa.tau) || length(prior.reg$learn.kappa.tau) != 1)
    stop("Argument learn.kappa.tau in regression prior must be a logical scalar.")
  if(is.null(prior.reg$learn.kappa.xi) || !is.logical(prior.reg$learn.kappa.xi) || length(prior.reg$learn.kappa.xi) != 1)
    stop("Argument learn.kappa.xi in regression prior must be a logical scalar.")
  if(!prior.reg$type %in% c("rw-t0", "rw-t1", "ind"))
    stop("Argument type in regression prior must either be 'rw-t0', 'rw-t1' or 'ind'.")
  check.hyper_double.positive(prior.reg$c, "c in regression prior")
  check.hyper_double.positive(prior.reg$B0, "B0 in regression prior")
  if(is.null(prior.reg$TG) || !is.logical(prior.reg$TG) || length(prior.reg$TG) != 1){
    stop("Argument TG in regression prior must be a logical scalar.")
  }
  if(is.null(prior.reg$TG.alternative) || !is.logical(prior.reg$TG.alternative) || length(prior.reg$TG.alternative) != 1){
    stop("Argument TG in regression prior must be a logical scalar.")
  }
  if(prior.reg$TG && prior.reg$type == "ind"){
    warning("Triple Gamma shrinkage will not be applied as you have selected an independence prior.")
  }
  if(prior.reg$TG && (prior.reg$a.tau >= 0.5 || prior.reg$a.tau <= 0 ||
                      prior.reg$a.xi >= 0.5 || prior.reg$a.xi <= 0 ||
                      prior.reg$c.tau >= 0.5 || prior.reg$c.tau <= 0 ||
                      prior.reg$c.xi >= 0.5 || prior.reg$c.xi <= 0)){
    warning("Arguments 'a.tau', 'a.xi', 'c.tau' and 'c.xi' have support on (0, 0.5) when using the triple Gamma prior and when learning those parameters! Setting them at fixed values, e.g., for Strawderman-Berger prior is fine.")
  }
  if(!prior.reg$TG && prior.reg$TG.alternative){
    stop("The alternative Triple Gamma priori requires argument TG to be TRUE as well.")
  }
}

check.prior.load <- function(prior.load){
  if(is.null(prior.load)) stop("Your loading prior setting is not allowed to be NULL.")
  check.hyper_double.positive(prior.load$d.phi, "d.phi in loading prior")
  check.hyper_double.positive(prior.load$e.phi, "e.phi in loading prior")
  check.hyper_double.positive(prior.load$d.zeta, "d.zeta in loading prior")
  check.hyper_double.positive(prior.load$e.zeta, "e.zeta in loading prior")
  check.hyper_double.positive(prior.load$a.phi, "a.phi in loading prior")
  check.hyper_double.positive(prior.load$kappa.phi, "kappa.phi in loading prior")
  check.hyper_double.positive(prior.load$a.zeta, "a.zeta in loading prior")
  check.hyper_double.positive(prior.load$kappa.zeta, "kappa.zeta in loading prior")
  if(is.null(prior.load$learn.kappa.phi) || !is.logical(prior.load$learn.kappa.phi) || length(prior.load$learn.kappa.phi) != 1)
    stop("Argument learn.kappa.phi in loading prior must be a logical scalar.")
  if(is.null(prior.load$learn.kappa.zeta) || !is.logical(prior.load$learn.kappa.zeta) || length(prior.load$learn.kappa.zeta) != 1)
    stop("Argument learn.kappa.zeta in loading prior must be a logical scalar.")
  if(is.null(prior.load$type) || !prior.load$type %in% c("rw-t0", "rw-t1", "ind", "cps"))
    stop("Argument type in in loading prior must either be 'rw-t0', 'rw-t1', 'ind' or 'cps'.")
  check.hyper_double.positive(prior.load$c, "c in loading prior")
  check.hyper_double.positive(prior.load$L0, "L0 in loading prior")
}

check.mcmc.opt <- function(mcmc.opt) {
  if(is.null(mcmc.opt)) stop("Your MCMC settings are not allowed to be NULL.")
  if (is.null(mcmc.opt$chain.length) || is.null(mcmc.opt$burnin) ||
      is.null(mcmc.opt$thin) || is.null(mcmc.opt$asis)) {
    stop("The arguments in 'mcmc.opt' are not allowed to be NULL.")
  }
  if (!is.numeric(mcmc.opt$chain.length) | length(mcmc.opt$chain.length) != 1) {
    stop("'mcmc.opt$chain.length' of wrong type or not single value")
  }
  if (!is.numeric(mcmc.opt$burnin) | length(mcmc.opt$burnin) != 1) {
    stop("'mcmc.opt$burn-in' of wrong type or not single value")
  }
  if (!is.logical(mcmc.opt$asis) | length(mcmc.opt$asis) != 1) {
    stop("'mcmc.opt$asis' of wrong type or not single value")
  }

  if (mcmc.opt$chain.length <= 0 | (floor(mcmc.opt$chain.length) != mcmc.opt$chain.length)) {
    stop("'mcmc.opt$chain.length' needs to be a positive integer")
  }
  if (mcmc.opt$burnin < 0 | (floor(mcmc.opt$burnin) != mcmc.opt$burnin)) {
    stop("'mcmc.opt$burn-in' needs to be a positive integer or zero")
  }
  if (mcmc.opt$burnin >= mcmc.opt$chain.length) {
    stop("burn-in period needs to be shorter than the entire chain length")
  }
  if (mcmc.opt$chain.length - mcmc.opt$burnin == 1){
    stop("Markov Chain needs to be at least of length 2 after burn-in")
  }
  if (mcmc.opt$thin > (mcmc.opt$chain.length - mcmc.opt$burnin)) {
    stop(paste("Thinning factor cannot be greater than the number of remaining samples after burn-in (",
               mcmc.opt$chain.length - mcmc.opt$burnin, ").", sep=""))
  }
  if ((is.numeric(mcmc.opt$thin) & !isTRUE(all.equal(mcmc.opt$thin, as.integer(mcmc.opt$thin)))) | mcmc.opt$thin <= 0) {
    stop("Thining factor must be a positive integer.")
  }
  if ((mcmc.opt$chain.length - mcmc.opt$burnin) %% mcmc.opt$thin != 0) {
    stop(paste("Thinning factor must be a multiple of remaining samples after burn-in ("),
         mcmc.opt$chain.length - mcmc.opt$burnin, ").", sep = "")
  }
}

check.settings.NegBin <- function(settings.NegBin){
  if(is.null(settings.NegBin)) stop("Your Negative Binomial settings are not allowed to be NULL for count models.")
  check.hyper_double.positive(settings.NegBin$alpha.r, "alpha.r in settings.NegBin")
  check.hyper_double.positive(settings.NegBin$beta.r, "beta.r in settings.NegBin")
  if(is.null(settings.NegBin$expansion.steps) || !is.numeric(settings.NegBin$expansion.steps) ||
     length(settings.NegBin$expansion.steps) != 1 || settings.NegBin$expansion.steps %% 1 != 0 ||
     !is.finite(settings.NegBin$expansion.steps) || settings.NegBin$expansion.steps < 1)
    stop("Argument 'expansion.steps' in settings.NegBin must be a single, positive integer.")
  check.hyper_double.positive(settings.NegBin$width, "width in settings.NegBin")
  if(is.null(settings.NegBin$p.overrelax) || !is.numeric(settings.NegBin$p.overrelax) ||
     settings.NegBin$p.overrelax < 0 || settings.NegBin$p.overrelax > 1)
    stop("Argument p.overrelax in settings.NegBin prior must be valid probability")
  if(settings.NegBin$p.overrelax > 0){
    check.hyper_double.positive(settings.NegBin$accuracy.overrelax, "accuracy.overrelax in settings.NegBin")
  }
}

# check input arguments of predict S3-functions --------------------------------

check.predict <- function(model, X.new, timepoint, coverage, pop.pred, n.replicates){

  if(!(is.matrix(X.new) || is.data.frame(X.new))){
    stop("Argument 'X.new' must be either a matrix or a data frame.")
  }

  if(ncol(X.new) != ncol(model$data$X)){
    stop("Number of columns in X.new must match number of columns in original design matrix. In case you fitted an intercept, X.new must contain a column of ones.")
  }

  if(is.null(colnames(X.new))){
    warning("Columns of argument 'X.new' are not labelled. This is dangerous as you must make sure that the columns match the columns of the original design matrix.")
  }

  else{

    if(any(!(colnames(X.new) %in% colnames(model$data$X)))){
      stop("Variable names of X.new must match the corresponding names from the original design matrix.")
    }

    if(sum(colnames(X.new) != colnames(model$data$X)) > 0){
      warning("Column names in 'X.new' were reordered to match the corresponding ones from the original design matrix.")
      X.new <- X.new[, colnames(model$data$X), drop = FALSE] # reordering variables w.r.t. original design matrix
    }

  }

  if(is.null(timepoint) || !is.numeric(timepoint) || length(timepoint) != 1 || timepoint %% 1 != 0 || !is.finite(timepoint)){
    stop("Argument 'timepoint' must be a single integer. If you want predictions for multiple time points, please call the predict function multiple times.")
  }

  if(!timepoint %in% unique(model$data$timeidx)){
    stop("Argument 'timepoint' must contain a time point for which the parameters were learned.")
  }

  if(is.null(coverage) || !is.numeric(coverage) || length(coverage) != 1 || !is.finite(coverage) || coverage < 0 || coverage > 1){
    stop("Argument 'coverage' must be a single numeric value between 0 and 1.")
  }

  if(is.null(pop.pred) || !is.logical(pop.pred) || length(pop.pred) != 1 || !is.finite(pop.pred)){
    stop("Argument 'pop.pred' must be a single logical value.")
  }

  if(is.null(n.replicates) || !is.numeric(n.replicates) || length(n.replicates) != 1 || timepoint %% 1 != 0 || !is.finite(n.replicates) || n.replicates < 1){
    stop("Argument 'n.replicates' must be a single, positive integer.")
  }

}

check.predict_ZINB <- function(model, X_nb.new, X_logit.new, timepoint, coverage, pop.pred, n.replicates){

  if(!(is.matrix(X_nb.new) || is.data.frame(X_nb.new))){
    stop("Argument 'X_nb.new' must be either a matrix or a data frame.")
  }

  if(ncol(X_nb.new) != ncol(model$data$X_nb)){
    stop("Number of columns in X_nb.new must match number of columns in original design matrix. In case you fitted an intercept, X_nb.new must contain a column of ones.")
  }

  if(is.null(colnames(X_nb.new))){
    warning("Columns of argument 'X_nb.new' are not labelled. This is dangerous as you must make sure that the columns match the columns of the original design matrix.")
  }

  else{

    if(any(!(colnames(X_nb.new) %in% colnames(model$data$X_nb)))){
      stop("Variable names of X_nb.new must match the corresponding names from the original design matrix.")
    }

    if(sum(colnames(X_nb.new) != colnames(model$data$X_nb)) > 0){
      warning("Column names in 'X_nb.new' were reordered to match the corresponding ones from the original design matrix.")
      X_nb.new <- X_nb.new[, colnames(model$data$X_nb), drop = FALSE] # reordering variables w.r.t. original design matrix
    }

  }

  if(!(is.matrix(X_logit.new) || is.data.frame(X_logit.new))){
    stop("Argument 'X_logit.new' must be either a matrix or a data frame.")
  }

  if(ncol(X_logit.new) != ncol(model$data$X_logit)){
    stop("Number of columns in X_logit.new must match number of columns in original design matrix. In case you fitted an intercept, X_logit.new must contain a column of ones.")
  }

  if(is.null(colnames(X_logit.new))){
    warning("Columns of argument 'X_logit.new' are not labelled. This is dangerous as you must make sure that the columns match the columns of the original design matrix.")
  }

  else{

    if(any(!(colnames(X_logit.new) %in% colnames(model$data$X_logit)))){
      stop("Variable names of X_logit.new must match the corresponding names from the original design matrix.")
    }

    if(sum(colnames(X_logit.new) != colnames(model$data$X_logit)) > 0){
      warning("Column names in 'X_logit.new' were reordered to match the corresponding ones from the original design matrix.")
      X_logit.new <- X_logit.new[, colnames(model$data$X_logit), drop = FALSE] # reordering variables w.r.t. original design matrix
    }

  }

  if(is.null(timepoint) || !is.numeric(timepoint) || length(timepoint) != 1 || timepoint %% 1 != 0 || !is.finite(timepoint)){
    stop("Argument 'timepoint' must be a single numeric value. If you want predictions for multiple time points, please call the predict function multiple times.")
  }

  if(is.null(timepoint) || !is.numeric(coverage) || length(coverage) != 1 || !is.finite(coverage) || coverage < 0 || coverage > 1){
    stop("Argument 'coverage' must be a single numeric value between 0 and 1.")
  }

  if(is.null(timepoint) || !is.logical(pop.pred) || length(pop.pred) != 1 || !is.finite(pop.pred)){
    stop("Argument 'pop.pred' must be a single logical value.")
  }

  if(is.null(timepoint) || !is.numeric(n.replicates) || length(n.replicates) != 1 || n.replicates %% 1 != 0 || !is.finite(n.replicates) || n.replicates < 1){
    stop("Argument 'n.replicates' must be a single, positive integer.")
  }

}
