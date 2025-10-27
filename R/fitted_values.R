compute_fitted_Gaussian_Probit_Logit_NegBin <- function(result){

  mcmc <- result$mcmc
  fmcmc <- result$fmcmc # M x n
  X <- result$data$X
  time <- result$data$timeidx
  Tmax <- max(time)
  id <- result$data$idx
  d <- ncol(X)
  N <- nrow(X)
  S <- nrow(mcmc)

  if(sum(startsWith(colnames(mcmc), "lambda_t")) == 1){ # cps
    la <- replicate(Tmax, mcmc[,"lambda_t"])
    colnames(la) <- paste0("lambda_t",1:Tmax)
    mcmc <- cbind(mcmc, la)
  }

  y.fit <- matrix(nrow = N, ncol = S)
  beta_list <- lapply(1:Tmax, function(t) {
    mcmc[, paste0("beta_t", 1:d, t), drop = FALSE]
  })
  lambda_list <- lapply(1:Tmax, function(t) {
    mcmc[, paste0("lambda_t", t)]
  })

  for(o in 1:N){
    t <- time[o]
    i <- id[o]
    x_row <- X[o, ]
    beta_t <- beta_list[[t]]
    f_i <- fmcmc[, i]
    lambda_t <- lambda_list[[t]]
    eta <- rowSums(beta_t * matrix(rep(x_row, each = S), nrow = S)) +
      f_i * lambda_t
    y.fit[o, ] <- sampling_fitted_auxiliary(model = result$model,
                                            mcmc = mcmc,
                                            eta = eta,
                                            S = S)
  }

  return(float::fl(y.fit))

}

compute_fitted_Gaussian_Probit_Logit_NegBin_no.fac <- function(result){

  mcmc <- result$mcmc
  X <- result$data$X
  time <- result$data$timeidx
  Tmax <- max(time)
  id <- result$data$idx
  d <- ncol(X)
  N <- nrow(X)
  S <- nrow(mcmc)

  y.fit <- matrix(nrow = N, ncol = S)
  beta_list <- lapply(1:Tmax, function(t) {
    mcmc[, paste0("beta_t", 1:d, t), drop = FALSE]
  })

  for(o in 1:N){
    t <- time[o]
    i <- id[o]
    x_row <- X[o, ]
    beta_t <- beta_list[[t]]
    eta <- rowSums(beta_t * matrix(rep(x_row, each = S), nrow = S))
    y.fit[o, ] <- sampling_fitted_auxiliary(model = result$model,
                                            mcmc = mcmc,
                                            eta = eta,
                                            S = S)
  }

  return(y.fit)

}

sampling_fitted_auxiliary <- function(model, mcmc, eta, S){
  if(model == "Gaussian"){
    sigma2 <- mcmc[, "sigma2"]
    y <- rnorm(S, mean = eta, sd = sqrt(sigma2))
  }
  else if(model == "Probit") y <- rbinom(S, size = 1, prob = pnorm(eta))
  else if(model == "Logit") y <- rbinom(S, size = 1, prob = plogis(eta))
  else if(model == "NegBin"){
    r <- mcmc[,"r"]
    y <- MASS::rnegbin(S, mu = r * exp(eta), theta = r)
  }
  return(y)
}

compute_fitted_ZINB <- function(result){

  # we do not need linear predictor of logit model as it is implicit in at-risk indicator w
  mcmc <- result$mcmc_nb
  fmcmc <- result$fmcmc_nb
  X <- result$data$X_nb
  time <- result$data$timeidx
  Tmax <- max(time)
  id <- result$data$idx
  d <- ncol(X)
  N <- nrow(X)
  S <- nrow(mcmc)

  if(sum(startsWith(colnames(mcmc), "lambda_t")) == 1){ # cps
    la <- replicate(Tmax, mcmc[,"lambda_t"])
    colnames(la) <- paste0("lambda_t",1:Tmax)
    mcmc <- cbind(mcmc, la)
  }

  y.fit <- matrix(nrow = N, ncol = S)
  beta_list <- lapply(1:Tmax, function(t) {
    mcmc[, paste0("beta_t", 1:d, t), drop = FALSE]
  })
  lambda_list <- lapply(1:Tmax, function(t) {
    mcmc[, paste0("lambda_t", t)]
  })

  r <- mcmc[,"r"]
  w <- result$mcmc_risk # S x N

  for(o in 1:N){
    t <- time[o]
    i <- id[o]
    x_row <- X[o, ]
    beta_t <- beta_list[[t]]
    lambda_t <- lambda_list[[t]]
    f_i <- fmcmc[,i]
    eta <- rowSums(beta_t * matrix(rep(x_row, each = S), nrow = S)) +
      f_i * lambda_t
    idx <- w[,o] == TRUE
    y.draw <- vector(mode = "numeric", length = length(w[,o]))
    y.draw[idx] <- MASS::rnegbin(sum(idx), mu = r[idx] * exp(eta[idx]), theta = r[idx])
    y.fit[o, ] <- y.draw
  }

  return(float::fl(y.fit))

}

compute_fitted_ZINB_no.fac <- function(result){

  # we do not need linear predictor of logit model as it is implicit in at-risk indicator w
  mcmc <- result$mcmc_nb
  X <- result$data$X_nb
  time <- result$data$timeidx
  Tmax <- max(time)
  id <- result$data$idx
  d <- ncol(X)
  N <- nrow(X)
  S <- nrow(mcmc)

  y.fit <- matrix(nrow = N, ncol = S)
  beta_list <- lapply(1:Tmax, function(t) {
    mcmc[, paste0("beta_t", 1:d, t), drop = FALSE]
  })

  r <- mcmc[,"r"]
  w <- result$mcmc_risk # S x N

  for(o in 1:N){
    t <- time[o]
    i <- id[o]
    x_row <- X[o, ]
    beta_t <- beta_list[[t]]
    eta <- rowSums(beta_t * matrix(rep(x_row, each = S), nrow = S))
    idx <- w[,o] == TRUE
    y.draw <- vector(mode = "numeric", length = length(w[,o]))
    y.draw[idx] <- MASS::rnegbin(sum(idx), mu = r[idx] * exp(eta[idx]), theta = r[idx])
    y.fit[o, ] <- y.draw
  }

  return(float::fl(y.fit))

}
