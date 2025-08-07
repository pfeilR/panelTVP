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

