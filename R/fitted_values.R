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

  return(y.fit)

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

  return(y.fit)

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

  return(y.fit)

}

compute_fitted_IV <- function(result){

  X1 <- result$data$X_stage1
  X2 <- result$data$X_stage2
  D_mcmc <- result$D_mcmc
  f_mcmc <- result$fmcmc_stage2
  mcmc1 <- result$mcmc_stage1
  mcmc2 <- result$mcmc_stage2
  rho_mcmc <- result$mcmc_rho
  N <- nrow(X1)
  S <- nrow(mcmc1)
  id <- result$data$idx
  time <- result$data$timeidx
  Tmax <- max(time)
  d1 <- ncol(X1)
  d2 <- ncol(X2)

  y.fit <- matrix(nrow = N, ncol = S)

  if(sum(startsWith(colnames(mcmc2), "lambda_t")) == 1){ # cps
    la <- replicate(Tmax, mcmc2[,"lambda_t"])
    colnames(la) <- paste0("lambda_t",1:Tmax)
    mcmc2 <- cbind(mcmc2, la)
  }

  beta1_list <- lapply(1:Tmax, function(t) {
    mcmc1[, paste0("beta_t", 1:d1, t), drop = FALSE]
  })
  beta2_list <- lapply(1:Tmax, function(t) {
    mcmc2[, paste0("beta_t", 1:d2, t), drop = FALSE]
  })
  lambda_list <- lapply(1:Tmax, function(t) {
    mcmc2[, paste0("lambda_t", t)]
  })

  for(o in 1:N){
    t <- time[o]
    i <- id[o]
    x1_row <- X1[o, ]
    x2_row <- X2[o,]
    beta1_t <- beta1_list[[t]]
    beta2_t <- beta2_list[[t]]
    lambda_t <- lambda_list[[t]]
    f_i <- f_mcmc[,i]

    eta1 <- rowSums(beta1_t * matrix(rep(x1_row, each = S), nrow = S))
    eta2 <- rowSums(beta2_t * matrix(rep(x2_row, each = S), nrow = S)) +
      f_i * lambda_t

    mu <- eta2 + rho_mcmc * sqrt(mcmc2[,"sigma2"]) * (D_mcmc[,o] - eta1)
    y.fit[o, ] <- rnorm(S, mean = mu, sd = sqrt(mcmc2[,"sigma2"] * (1 - rho_mcmc^2)))

  }

  return(y.fit)

}

compute_fitted_IV_no.fac <- function(result){

  X1 <- result$data$X_stage1
  X2 <- result$data$X_stage2
  D_mcmc <- result$D_mcmc
  mcmc1 <- result$mcmc_stage1
  mcmc2 <- result$mcmc_stage2
  rho_mcmc <- result$mcmc_rho
  N <- nrow(X1)
  S <- nrow(mcmc1)
  id <- result$data$idx
  time <- result$data$timeidx
  Tmax <- max(time)
  d1 <- ncol(X1)
  d2 <- ncol(X2)

  y.fit <- matrix(nrow = N, ncol = S)

  beta1_list <- lapply(1:Tmax, function(t) {
    mcmc1[, paste0("beta_t", 1:d1, t), drop = FALSE]
  })
  beta2_list <- lapply(1:Tmax, function(t) {
    mcmc2[, paste0("beta_t", 1:d2, t), drop = FALSE]
  })

  for(o in 1:N){
    t <- time[o]
    i <- id[o]
    x1_row <- X1[o, ]
    x2_row <- X2[o,]
    beta1_t <- beta1_list[[t]]
    beta2_t <- beta2_list[[t]]

    eta1 <- rowSums(beta1_t * matrix(rep(x1_row, each = S), nrow = S))
    eta2 <- rowSums(beta2_t * matrix(rep(x2_row, each = S), nrow = S))

    mu <- eta2 + rho_mcmc * sqrt(mcmc2[,"sigma2"]) * (D_mcmc[,o] - eta1)
    y.fit[o, ] <- rnorm(S, mean = mu, sd = sqrt(mcmc2[,"sigma2"] * (1 - rho_mcmc^2)))

  }

  return(y.fit)

}

