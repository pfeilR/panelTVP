sim_panelTVP_IV <- function(n,
                            Tmax,
                            beta_stage1,
                            theta_stage1,
                            lambda_stage1,
                            psi_stage1,
                            beta_stage2,
                            theta_stage2,
                            lambda_stage2,
                            beta_D,
                            theta_D,
                            psi_stage2,
                            rho,
                            sigma2,
                            n.instruments = 1){

  # Input Checks ---------------------------------------------------------------

  # NEEDED HERE !!!

  t <- Tmax
  d1 <- length(beta_stage1)
  d2 <- length(beta_stage2)

  # STAGE 1 --------------------------------------------------------------------

  # beta
  beta0_stage1 <- rnorm(d1)
  betat_nc_stage1 <- apply(matrix(c(beta0_stage1, rnorm(t*d1)), byrow = TRUE, ncol = d1), 2, cumsum)
  Dmat1 <- diag(theta_stage1)
  beta_t_stage1 <- matrix(beta_stage1, nrow = t+1, ncol = d1, byrow = TRUE) +
    t(apply(betat_nc_stage1, 1, function(x) Dmat1 %*% x))
  beta_t_stage1 <- beta_t_stage1[-1, , drop = FALSE]

  # lambda
  lambda_t_stage1 <- matrix(lambda_stage1 + psi_stage1 * rnorm(t), ncol = 1)

  # design matrices
  X_stage1 <- matrix(rnorm(n*t*d1), n*t, d1)
  X_stage1 <- model.matrix(~X_stage1[,-1])
  linpred_stage1 <- numeric(n*t)
  F_stage1 <- numeric(n*t)
  fi_stage1 <- base::scale(rnorm(n))
  for(tt in 1:t){
    ind <- ((tt-1)*n + 1):(n*tt)
    linpred_stage1[ind] <- X_stage1[ind, ] %*% beta_t_stage1[tt, ]
    F_stage1[ind] <- fi_stage1 * lambda_t_stage1[tt]
  }

  # STAGE 2 --------------------------------------------------------------------

  # beta
  beta0_stage2 <- rnorm(d2)
  betat_nc_stage2 <- apply(matrix(c(beta0_stage2, rnorm(t*d2)), byrow = TRUE, ncol = d2), 2, cumsum)
  Dmat2 <- diag(theta_stage2)
  beta_t_stage2 <- matrix(beta_stage2, nrow = t+1, ncol = d2, byrow = TRUE) +
    t(apply(betat_nc_stage2, 1, function(x) Dmat2 %*% x))
  beta_t_stage2 <- beta_t_stage2[-1, , drop = FALSE]

  # beta for treatment variable D
  beta0_D <- rnorm(1)
  betat_nc_D <- apply(matrix(c(beta0_D, rnorm(t)), byrow = TRUE, ncol = 1), 2, cumsum)
  beta_t_D <- matrix(matrix(beta_D, nrow = t+1, ncol = 1, byrow = TRUE) + theta_D * betat_nc_D, ncol = 1)
  beta_t_D <- beta_t_D[-1, , drop = FALSE]

  # lambda
  lambda_t_stage2 <- matrix(lambda_stage2 + psi_stage2 * rnorm(t), ncol = 1)

  # design matrices
  X_stage2 <- matrix(rnorm(n*t*d2), n*t, d2)
  X_stage2 <- model.matrix(~X_stage2[,-1])
  linpred_stage2 <- numeric(n*t)
  F_stage2 <- numeric(n*t)
  fi_stage2 <- base::scale(rnorm(n))
  for(tt in 1:t){
    ind <- ((tt-1)*n + 1):(n*tt)
    linpred_stage2[ind] <- X_stage2[ind, ] %*% beta_t_stage2[tt, ]
    F_stage2[ind] <- fi_stage2 * lambda_t_stage2[tt]
  }

  # Error distribution
  Sigma <- matrix(c(1, rho*sqrt(sigma2),
                    rho*sqrt(sigma2), sigma2), nrow = 2, byrow = TRUE)
  eps <- MASS::mvrnorm(n*t, mu = c(0,0), Sigma = Sigma)
  epsilon1 <- eps[,1]
  epsilon2 <- eps[,2]

  # Treatment simulation -------------------------------------------------------

  D_star <- linpred_stage1 + F_stage1 + epsilon1
  D <- ifelse(D_star > 0, 1, 0)
  D_ <- numeric(n*t)
  for(tt in 1:t){
    ind <- ((tt-1)*n + 1):(n*tt)
    D_[ind] <- D[ind] * beta_t_D[tt,]
  }

  # Outcome simulation ---------------------------------------------------------

  y <- linpred_stage2 + F_stage2 + D_ + epsilon2

  # return data and parameters -------------------------------------------------

  id <- rep(1:n, times = t)
  t <- rep(1:t, each = n)

  X_stage1 <- as.matrix(X_stage1[,-1])
  # if((ncol(X_stage1)-1)==1) X_stage1 <- as.matrix(X_stage1)
  if(ncol(X_stage1) == 1) colnames(X_stage1) <- "X_stage1.Z1"
  else colnames(X_stage1) <- c(paste0("W", 1:(d1-n.instruments-1)), paste0("Z", 1:n.instruments))

  X_stage2 <- X_stage2[,-1]
  if((d2-1)==1) X_stage2 <- as.matrix(X_stage2)
  colnames(X_stage2) <- paste0("W", 1:(d2-1))

  observed <- data.frame(y = y, X_stage1 = X_stage1, X_stage2 = X_stage2, D = D,
                         t = t, id = id)
  ret <- list(observed = observed,
              D_star = D_star,
              treatment_effect = beta_t_D,
              beta_t_stage1 = beta_t_stage1,
              lambda_t_stage1 = lambda_t_stage1,
              beta_t_stage2 = beta_t_stage2,
              lambda_t_stage2 = lambda_t_stage2,
              epsilon = data.frame(epsilon1 = epsilon1, epsilon2 = epsilon2),
              sigma2 = sigma2,
              rho = rho)

  return(ret)

}
