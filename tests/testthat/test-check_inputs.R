# checks for formula -----------------------------------------------------------

test_that("missing formula", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "Gaussian"),
               "Arguments 'formula' and 'data' are required with no defaults.")
})

test_that("not a formula object", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(f = "Time to Wake up Little Girl",
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "Gaussian"),
               "Argument 'formula' must be a formula object.")
})

test_that("wrong formula for given model", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(y ~ W1 | W1,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "Gaussian"),
               "The character '|' in argument 'formula' is only valid for Zero-Inflated Negative Binomial regression.")
})

test_that("wrong formula for given model", {
  d <- sim_panelTVP(n = 1000,
                    Tmax = 4,
                    beta.nb = c(0.05,0),
                    theta.nb = c(0,0),
                    lambda.nb = 0.3,
                    psi.nb = 0,
                    beta.logit = c(0.01,0),
                    theta.logit = c(0.05,0),
                    lambda.logit = 0.2,
                    psi.logit = 0.05,
                    r = 0.5,
                    model = "ZINB")
  expect_error(panelTVP(y ~ W1.nb,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "ZINB"),
               "Argument 'formula' needs exactly one '|' character for Zero-Inflated Negative Binomial regression.")
})

test_that("formula contains variables not included in the dataset", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(y ~ W4,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "Gaussian"),
               "There are variables in your 'formula' argument that are not contained in 'data'.")
})

test_that("response also in linear predictor", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(W1 ~ W1,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "Gaussian"),
               "Response variable is also contained as predictor.")
})

test_that("response also in linear predictor (ZINB)", {
  d <- sim_panelTVP(n = 1000,
                    Tmax = 4,
                    beta.nb = c(0.05,0),
                    theta.nb = c(0,0),
                    lambda.nb = 0.3,
                    psi.nb = 0,
                    beta.logit = c(0.01,0),
                    theta.logit = c(0.05,0),
                    lambda.logit = 0.2,
                    psi.logit = 0.05,
                    r = 0.5,
                    model = "ZINB")
  expect_error(panelTVP(y ~ y | W1.logit,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "ZINB"),
               "Response variable is also contained as predictor.")
})

# model and response type ------------------------------------------------------

test_that("invalid model type", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(y ~ W1,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "Darth Maul"),
               "Argument 'model' must be either 'Gaussian', 'Probit', 'Logit', 'NegBin' or 'ZINB'.")
})

test_that("wrong response for binary model", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(y ~ W1,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "Logit"),
               "When response is of type numeric, it must only contain 0 and 1 for Probit and Logit models.")
})

test_that("wrong response for count model", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(y ~ W1,
                        data = d$observed,
                        id = d$observed$id,
                        t = d$observed$t,
                        model = "NegBin"),
               "Response must be positive for count data regression.")
})

# id and t ---------------------------------------------------------------

test_that("argument t not specified", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(y ~ W1,
                        data = d$observed,
                        id = d$observed$id,
                        model = "Gaussian"),
               "Argument 't' must be an integer-valued vector with length equal to the number of observations in 'data'.")
})

test_that("argument t of wrong dimension", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  d$observed$t <- NULL
  expect_error(panelTVP(y ~ W1,
                        data = d$observed,
                        id = d$observed$id,
                        t = 1:10,
                        model = "Gaussian"),
               "Argument 't' must be an integer-valued vector with length equal to the number of observations in 'data'.")
})

# -> es wird kein Fehler geworfen, warum nicht ???



test_that("argument id not specified", {
  d <- sim_panelTVP(n = 5000,
                    Tmax = 4,
                    beta = c(-1,1),
                    theta = c(0,0.5),
                    lambda = 2,
                    psi = 0.5,
                    model = "Gaussian",
                    sigma2 = 0.5)
  expect_error(panelTVP(y ~ W1,
                        data = d$observed,
                        t = d$observed$t,
                        model = "Gaussian"),
               "Argument 'id' must be an integer-valued vector with length equal to the number of observations in 'data'.")
})



check.panelTVP <- function(formula, data, id, t, model, prior.reg, prior.var, prior.load,
                           prior.reg_nb, prior.load_nb, prior.reg_logit, prior.load_logit,
                           mcmc.opt, settings.NegBin, HPD.coverage, R.WAIC,
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
  vars <- all.vars(formula)
  if(sum(!(vars %in% colnames(data))) > 0){
    stop("There are variables in your 'formula' argument that are not contained in 'data'.")
  }
  if(is.null(t) || length(t) != nrow(data) || !is.numeric(t) || sum(t %% 1) != 0 || sum(!is.finite(t)) != 0){
    stop("Argument 't' must be an integer-valued vector with length equal to the number of observations in 'data'.")
  }
  if(is.null(id) || length(id) != nrow(data) || !is.numeric(id) || sum(id %% 1) != 0 || sum(!is.finite(id)) != 0){
    stop("Argument 'id' must be an integer-valued vector with length equal to the number of observations in 'data'.")
  }

  # check if panel is balanced
  if(length(t) %% max(t) != 0){
    stop("Panel must be balanced, i.e., the same number of repeated measurements for every subject.")
  }

  # response variable check
  resp <- data[, as.character(formula[[2]])]
  if(model %in% c("Probit", "Logit")){
    if(is.factor(resp) && length(levels(resp)) != 2){
      stop("When response variable is a factor, it must have exactly two levels for Probit and Logit models.")
    } else if(is.numeric(resp) && !all(resp %in% c(0,1))){
      stop("When response is of type numeric, it must only contain 0 and 1 for Probit and Logit models.")
    } else if(!is.factor(resp) && !is.numeric(resp) && !is.logical(resp)){
      stop("Response must be a factor with two levels, numeric 0/1, or logical for Probit and Logit models.")
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

    check.prior.reg(prior.reg)
    if(model == "Gaussian"){
      # prior.var.check
      check.hyper_double.positive(prior.var$learn.C0.hyp$g0, "g0 in prior.var")
      check.hyper_double.positive(prior.var$learn.C0.hyp$G0, "G0 in prior.var")
      check.hyper_double.positive(prior.var$c0, "c0 in prior.var")
    }
    check.prior.load(prior.load)

  } else{ # ZINB-Checks

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
     !is.finite(R.WAIC) || R.WAIC < 1)
    stop("Argument 'R.WAIC' must be a single, positive integer.")

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
