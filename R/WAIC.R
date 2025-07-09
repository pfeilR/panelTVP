# In this file, the WAIC is computed under 6 different likelihoods.

library(LaplacesDemon)

# compute linear predictor based on model output
# this will return a matrix of dimension (nT x S),
# where S = (chain.length - burnin)/thin, i.e., the remaining draws of model output

# Update: 12.06.2025 (code adapted to allow for cps priori)
# Last update: 24.06.25 (omitting all rows with NAs in response)
# this is important as WAIC is based soley on likelihoods of observables!!!
# (see also how brms handles this!)

lp.model <- function(model){

  X <- cbind(model$data$X, t = model$data$timeidx)
  fi <- model$fmcmc
  d <- model$data$d
  Tmax <- model$data$Tmax
  mc <- model$mcmc
  S <- nrow(mc)
  N <- nrow(X)
  if(sum(startsWith(colnames(mc), "lambda_t")) == 1) cps <- TRUE
  else cps <- FALSE
  eta <- matrix(NA, nrow = N, ncol = S)
  X_split <- split(as.data.frame(X), X[,"t"])
  for(s in 1:S){
    eta_s <- vector("list", Tmax)
    if(cps) lambda_col <- "lambda_t"
    for(t in 1:Tmax){
      Xt <- as.matrix(X_split[[as.character(t)]][, -ncol(X)])
      beta_cols <- paste0("beta_t", 1:d, t)
      if(!cps) lambda_col <- paste0("lambda_t", t)
      eta_s[[t]] <- Xt %*% mc[s, beta_cols] + fi[s,] * mc[s, lambda_col]
    }
    eta[,s] <- do.call("rbind", eta_s)
  }

  return(eta)

}

lp.zinb <- function(model){

  Tmax <- model$data$Tmax

  # logit
  X_logit <- cbind(model$data$X_logit, t = model$data$timeidx)
  fi_logit <- model$fmcmc_logit
  d_logit <- model$data$d_logit
  mc_logit <- model$mcmc_logit
  S <- nrow(mc_logit)
  N <- nrow(X_logit)
  if(sum(startsWith(colnames(mc_logit), "lambda_t")) == 1) cps_logit <- TRUE
  else cps_logit <- FALSE
  if(cps_logit) lambda_col_logit <- "lambda_t"
  eta_logit <- matrix(NA, nrow = N, ncol = S)
  X_split_logit <- split(as.data.frame(X_logit), X_logit[,"t"])
  for(s in 1:S){
    eta_s_logit <- vector("list", Tmax)
    for(t in 1:Tmax){
      Xt_logit <- as.matrix(X_split_logit[[as.character(t)]][, -ncol(X_logit)])
      beta_cols_logit <- paste0("beta_t", 1:d_logit, t)
      if(!cps_logit) lambda_col_logit <- paste0("lambda_t", t)
      eta_s_logit[[t]] <- Xt_logit %*% mc_logit[s, beta_cols_logit] +
        fi_logit[s,] * mc_logit[s, lambda_col_logit]
    }
    eta_logit[,s] <- do.call("rbind", eta_s_logit)
  }

  # negative binomial
  X_nb <- cbind(model$data$X_nb, t = model$data$timeidx)
  fi_nb <- model$fmcmc_nb
  d_nb <- model$data$d_nb
  mc_nb <- model$mcmc_nb
  S <- nrow(mc_nb)
  N <- nrow(X_nb)
  if(sum(startsWith(colnames(mc_nb), "lambda_t")) == 1) cps_nb <- TRUE
  else cps_nb <- FALSE
  if(cps_nb) lambda_col_nb <- "lambda_t"
  eta_nb <- matrix(NA, nrow = N, ncol = S)
  X_split_nb <- split(as.data.frame(X_nb), X_nb[,"t"])
  for(s in 1:S){
    eta_s_nb <- vector("list", Tmax)
    for(t in 1:Tmax){
      Xt_nb <- as.matrix(X_split_nb[[as.character(t)]][, -ncol(X_nb)])
      beta_cols_nb <- paste0("beta_t", 1:d_nb, t)
      if(!cps_nb) lambda_col_nb <- paste0("lambda_t", t)
      eta_s_nb[[t]] <- Xt_nb %*% mc_nb[s, beta_cols_nb] +
        fi_nb[s,] * mc_nb[s, lambda_col_nb]
    }
    eta_nb[,s] <- do.call("rbind", eta_s_nb)
  }

  etas <- list(eta_logit = eta_logit, eta_nb = eta_nb)

  return(etas)

}

compute_waic <- function(model){

  m <- model$model
  idx.observed <- !is.na(model$data$y)
  y <- model$data$y[idx.observed]

  if(m != "ZINB"){

  eta <- lp.model(model)[idx.observed,]
  S <- ncol(eta)
  ll <- matrix(NA, nrow = length(y), ncol = S)

  if(m == "Gaussian"){
    sigma2 <- model$mcmc[,"sigma2"]
    for(s in 1:S){
      ll[,s] <- dnorm(y, mean = eta[,s], sd = sqrt(sigma2[s]), log = TRUE)
    }
  }

  if(m == "Probit"){
    p <- pnorm(eta)
    p <- pmax(pmin(p, 1 - 1e-4), 1e-4)
    for(s in 1:S){
      ll[,s] <- dbinom(y, size = 1, prob = p[,s], log = TRUE)
    }
  }

  if(m == "Logit"){
    p <- plogis(eta)
    p <- pmax(pmin(p, 1 - 1e-4), 1e-4)
    for(s in 1:S){
      ll[,s] <- dbinom(y, size = 1, prob = p[,s], log = TRUE)
    }
  }

  if(m == "NegBin"){
    r <- model$mcmc[,"r"]
    for(s in 1:S){
      ll[,s] <- dnbinom(y, size = r[s], mu = r[s] * exp(eta[,s]), log = TRUE)
    }
  }

  } else{ # ZINB model

    etas <- lp.zinb(model)
    eta_logit <- etas[["eta_logit"]][idx.observed,]
    eta_nb <- etas[["eta_nb"]][idx.observed,]
    S <- ncol(eta_logit)
    p.risk <- pmax(pmin(plogis(eta_logit), 1 - 1e-4), 1e-4)
    r <- model$mcmc_nb[,"r"]
    yS <- matrix(rep(y, S), ncol = S)
    rS <- matrix(r, nrow = length(y), ncol = S, byrow = TRUE)
    nb.lik <- dnbinom(yS, size = rS, mu = rS * exp(eta_nb))
    ll <- matrix(NA, nrow = length(y), ncol = S)
    zeros <- y == 0
    for(s in 1:S){
      w <- model$mcmc_risk[s,idx.observed]
      # structural zero
      ll[zeros & w == 0, s] <- log(1 - p.risk[zeros & w == 0, s])
      # at-risk zero
      ll[zeros & w == 1, s] <- log(p.risk[zeros & w == 1, s] * nb.lik[zeros & w == 1, s])
      # at risk with positive count
      ll[!zeros & w == 1, s] <- log(p.risk[!zeros & w == 1, s] * nb.lik[!zeros & w == 1, s])
    }

  }

  waic <- LaplacesDemon::WAIC(ll)$WAIC

  if(sum(is.na(model$data$y)) > 0){
    warning("NAs are present in response variable. WAIC was computed based on the observed data.")
  }

  return(waic)

}
