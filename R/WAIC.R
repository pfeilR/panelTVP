# In this file, the WAIC is computed under 5 different likelihoods.

# compute linear predictor based on model output
# this will return a matrix of dimension (nT x S),
# where S = (chain.length - burnin)/thin, i.e., the remaining draws of model output

# Update: 12.06.2025 (code adapted to allow for cps priori)
# Last update: 24.06.25 (omitting all rows with NAs in response)
# this is important as WAIC is based soley on likelihoods of observables!!!
# (see also how brms handles this!)

# Main function ----------------------------------------------------------------

compute_waic <- function(model, random.effects, R){

  m <- model$model
  idx.observed <- !is.na(model$data$y)
  y <- model$data$y[idx.observed]

  if(m != "ZINB"){

    S <- nrow(model$mcmc)
    ll <- matrix(NA, nrow = length(y), ncol = S)

    if(!random.effects){
      eta <- lp.model_no.fac(model)[idx.observed,]
      for(s in 1:S){
        ll[,s] <- log_lik_army(y, eta[,s], model, s)
      }
    } else{
      for(s in 1:S){
        eta <- lp.model(model, s = s, R = R)[idx.observed,]
        l_r <- apply(eta, 2, function(x){
          exp(log_lik_army(y, x, model, s))
        })
        ll[,s] <- log(rowMeans(l_r))
      }
    }

  } else{ # ZINB model

    S <- nrow(model$mcmc_logit)
    r <- model$mcmc_nb[,"r"]
    ll <- matrix(NA, nrow = length(y), ncol = S)

    if(!random.effects){
      etas <- lp.zinb_no.fac(model)
      eta_logit <- etas[["eta_logit"]][idx.observed,]
      eta_nb <- etas[["eta_nb"]][idx.observed,]
      p.risk <- apply(eta_logit, 2, function(x) pmax(pmin(plogis(x), 1 - 1e-4), 1e-4))
      for(s in 1:S){
        ll[,s] <- log(ifelse(y > 0,
                         p.risk[,s] * dnbinom(y, size = r[s], mu = r[s] * exp(eta_nb[,s])),
                         (1 - p.risk[,s]) + p.risk[,s] * dnbinom(0, size = r[s], mu = r[s] * exp(eta_nb[,s]))
                  ))
      }

    } else{ # random.effects = TRUE

      for(s in 1:S){
        etas <- lp.zinb(model, s = s, R = R)
        eta_logit <- etas[["eta_logit"]][idx.observed,]
        eta_nb <- etas[["eta_nb"]][idx.observed,]
        p.risk <- apply(eta_logit, 2, function(x) pmax(pmin(plogis(x), 1 - 1e-4), 1e-4))
        l_r <- sapply(1:R, function(replicant){
          mu <- r[s] * exp(eta_nb[,replicant])
          pr <- p.risk[,replicant]
          ifelse(y > 0,
                 pr * dnbinom(y, size = r[s], mu = mu),
                 (1 - pr) + pr * dnbinom(0, size = r[s], mu = mu))
        })
        ll[,s] <- log(rowMeans(l_r))
      }

    }

  }

  waic <- LaplacesDemon::WAIC(ll)$WAIC

  if(sum(is.na(model$data$y)) > 0){
    warning("NAs are present in response variable. WAIC was computed based on the observed data.")
  }

  return(waic)

}

# Helper functions -------------------------------------------------------------

log_lik_army <- function(y, eta, model, s){

  m <- model$model

  if(m == "Gaussian"){
    sigma2 <- model$mcmc[s, "sigma2"]
    return(dnorm(y, mean = eta, sd = sqrt(sigma2), log = TRUE))
  } else if(m == "Probit"){
    p <- pnorm(eta)
    p <- pmax(pmin(p, 1 - 1e-4), 1e-4)
    return(dbinom(y, size = 1, prob = p, log = TRUE))
  } else if(m == "Logit"){
    p <- plogis(eta)
    p <- pmax(pmin(p, 1 - 1e-4), 1e-4)
    return(dbinom(y, size = 1, prob = p, log = TRUE))
  } else if (m == "NegBin"){
    r <- model$mcmc[s, "r"]
    mu <- r * exp(eta)
    return(dnbinom(y, size = r, mu = mu, log = TRUE))
  } else{
    stop("Not supported model type.")
  }

}

lp.model <- function(model, s, R){

  X <- cbind(model$data$X, t = model$data$timeidx)
  d <- model$data$d
  Tmax <- model$data$Tmax
  mc <- model$mcmc
  N <- nrow(X)
  if(sum(startsWith(colnames(mc), "lambda_t")) == 1) cps <- TRUE
  else cps <- FALSE
  X_split <- split(as.data.frame(X), X[,"t"])
  if(cps) lambda_col <- "lambda_t"
  eta_list <- vector("list", Tmax)
  fi_matrix <- matrix(rnorm((N/Tmax)*R), ncol = R)
  fi_matrix <- scale(fi_matrix)
  for(t in 1:Tmax){
    Xt <- as.matrix(X_split[[as.character(t)]][, -ncol(X)])
    beta_cols <- paste0("beta_t", 1:d, t)
    if(!cps) lambda_col <- paste0("lambda_t", t)
    lambda <- as.numeric(mc[s, lambda_col])
    eta_fix <- Xt %*% mc[s, beta_cols]
    eta_mat <- replicate(R, c(eta_fix)) + fi_matrix * lambda
    eta_list[[t]] <- eta_mat
  }

  eta <- do.call(rbind, eta_list)  # gives me an (N × R) matrix
  return(eta)

}

lp.model_no.fac <- function(model){

  X <- cbind(model$data$X, t = model$data$timeidx)
  d <- model$data$d
  Tmax <- model$data$Tmax
  mc <- model$mcmc
  S <- nrow(mc)
  N <- nrow(X)
  eta <- matrix(NA, nrow = N, ncol = S)
  X_split <- split(as.data.frame(X), X[,"t"])
  for(s in 1:S){
    eta_s <- vector("list", Tmax)
    for(t in 1:Tmax){
      Xt <- as.matrix(X_split[[as.character(t)]][, -ncol(X)])
      beta_cols <- paste0("beta_t", 1:d, t)
      eta_s[[t]] <- Xt %*% mc[s, beta_cols]
    }
    eta[,s] <- do.call("rbind", eta_s)
  }

  return(eta)

}

lp.zinb <- function(model, s, R){

  X_logit <- cbind(model$data$X_logit, t = model$data$timeidx)
  d_logit <- model$data$d_logit
  Tmax <- model$data$Tmax
  mc_logit <- model$mcmc_logit
  N <- nrow(X_logit)
  if(sum(startsWith(colnames(mc_logit), "lambda_t")) == 1) cps_logit <- TRUE
  else cps_logit <- FALSE
  X_split_logit <- split(as.data.frame(X_logit), X_logit[,"t"])
  if(cps_logit) lambda_col_logit <- "lambda_t"
  eta_list_logit <- vector("list", Tmax)
  fi_matrix_logit <- matrix(rnorm((N/Tmax)*R), ncol = R)
  fi_matrix_logit <- scale(fi_matrix_logit)

  X_nb <- cbind(model$data$X_nb, t = model$data$timeidx)
  d_nb <- model$data$d_nb
  Tmax <- model$data$Tmax
  mc_nb <- model$mcmc_nb
  N <- nrow(X_nb)
  if(sum(startsWith(colnames(mc_nb), "lambda_t")) == 1) cps_nb <- TRUE
  else cps_nb <- FALSE
  X_split_nb <- split(as.data.frame(X_nb), X_nb[,"t"])
  if(cps_nb) lambda_col_nb <- "lambda_t"
  eta_list_nb <- vector("list", Tmax)
  fi_matrix_nb <- matrix(rnorm((N/Tmax)*R), ncol = R)
  fi_matrix_nb <- scale(fi_matrix_nb)

  for(t in 1:Tmax){

    Xt_logit <- as.matrix(X_split_logit[[as.character(t)]][, -ncol(X_logit)])
    beta_cols_logit <- paste0("beta_t", 1:d_logit, t)
    if(!cps_logit) lambda_col_logit <- paste0("lambda_t", t)
    lambda_logit <- as.numeric(mc_logit[s, lambda_col_logit])
    eta_fix_logit <- Xt_logit %*% mc_logit[s, beta_cols_logit]
    eta_mat_logit <- replicate(R, c(eta_fix_logit)) + fi_matrix_logit * lambda_logit
    eta_list_logit[[t]] <- eta_mat_logit

    Xt_nb <- as.matrix(X_split_nb[[as.character(t)]][, -ncol(X_nb)])
    beta_cols_nb <- paste0("beta_t", 1:d_nb, t)
    if(!cps_nb) lambda_col_nb <- paste0("lambda_t", t)
    lambda_nb <- as.numeric(mc_nb[s, lambda_col_nb])
    eta_fix_nb <- Xt_nb %*% mc_nb[s, beta_cols_nb]
    eta_mat_nb <- replicate(R, c(eta_fix_nb)) + fi_matrix_nb * lambda_nb
    eta_list_nb[[t]] <- eta_mat_nb

  }

  eta_logit <- do.call(rbind, eta_list_logit)
  eta_nb <- do.call(rbind, eta_list_nb)  # gives me an (N × R) matrix

  etas <- list(eta_logit = eta_logit, eta_nb = eta_nb)
  return(etas)

}

lp.zinb_no.fac <- function(model){

  Tmax <- model$data$Tmax

  # logit
  X_logit <- cbind(model$data$X_logit, t = model$data$timeidx)
  d_logit <- model$data$d_logit
  mc_logit <- model$mcmc_logit
  S <- nrow(mc_logit)
  N <- nrow(X_logit)
  eta_logit <- matrix(NA, nrow = N, ncol = S)
  X_split_logit <- split(as.data.frame(X_logit), X_logit[,"t"])
  for(s in 1:S){
    eta_s_logit <- vector("list", Tmax)
    for(t in 1:Tmax){
      Xt_logit <- as.matrix(X_split_logit[[as.character(t)]][, -ncol(X_logit)])
      beta_cols_logit <- paste0("beta_t", 1:d_logit, t)
      eta_s_logit[[t]] <- Xt_logit %*% mc_logit[s, beta_cols_logit]
    }
    eta_logit[,s] <- do.call("rbind", eta_s_logit)
  }

  # negative binomial
  X_nb <- cbind(model$data$X_nb, t = model$data$timeidx)
  d_nb <- model$data$d_nb
  mc_nb <- model$mcmc_nb
  S <- nrow(mc_nb)
  N <- nrow(X_nb)
  eta_nb <- matrix(NA, nrow = N, ncol = S)
  X_split_nb <- split(as.data.frame(X_nb), X_nb[,"t"])
  for(s in 1:S){
    eta_s_nb <- vector("list", Tmax)
    for(t in 1:Tmax){
      Xt_nb <- as.matrix(X_split_nb[[as.character(t)]][, -ncol(X_nb)])
      beta_cols_nb <- paste0("beta_t", 1:d_nb, t)
      eta_s_nb[[t]] <- Xt_nb %*% mc_nb[s, beta_cols_nb]
    }
    eta_nb[,s] <- do.call("rbind", eta_s_nb)
  }

  etas <- list(eta_logit = eta_logit, eta_nb = eta_nb)

  return(etas)

}
