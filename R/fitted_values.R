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
  li <- fitted_identify_lambda(lambda_list)
  lambda_list <- li$identified_lambda
  sign.lambda1 <- li$sign.lambda1

  for(o in 1:N){
    t <- time[o]
    i <- id[o]
    x_row <- X[o, ]
    beta_t <- beta_list[[t]]
    f_i <- fmcmc[, i]
    k <- kmeans(f_i, centers = 2)
    if(which.max(k$centers) == 1) f_i <- ifelse(f_i < 0, f_i, -f_i)
    else f_i <- ifelse(f_i > 0, f_i, -f_i)
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

fitted_identify_lambda <- function(lambda_list){

  sign.lambda1 <- sign(lambda_list[[1]])
  identified_lambda <- lapply(lambda_list, function(x) sign.lambda1*x)
  return(list(identified_lambda = identified_lambda, sign.lambda1 = sign.lambda1))

}
