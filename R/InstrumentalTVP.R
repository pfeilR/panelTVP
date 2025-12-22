InstrumentalTVP <- function(df,
                            prior.reg_stage1,
                            prior.reg_stage2,
                            prior.var_stage2,
                            prior.load_stage2,
                            prior.rho,
                            mcmc.opt,
                            alpha_stage1,
                            alpha_stage2,
                            lambda_stage2,
                            alpha_lambda_stage2,
                            reff_stage2,
                            tv.load_stage2,
                            res_frame_stage1,
                            res_frame_stage2,
                            f_sum_stage2,
                            f_mat_stage2,
                            miss,
                            sigma2v,
                            rho,
                            HPD.coverage,
                            random.effects,
                            progress.bar,
                            Treatment.Variable,
                            latent.utility.targeted){

  # return object for saving of rho
  rho_vec <- numeric(mcmc.opt$chain.length)

  # return object for saving D.star (only saved for fitted values)
  D_save <- matrix(nrow = mcmc.opt$chain.length, ncol = length(df$y))

  fi.count_stage2 <- 1

  df.stage2 <- df
  names(df.stage2)[names(df.stage2) == "X_stage2"] <- "X"
  names(df.stage2)[names(df.stage2) == "d_stage2"] <- "d"
  X.t_stage2 <- cbind(df.stage2$X, t = df.stage2$timeidx)
  df.stage2[["X_stage1"]] <- NULL
  df.stage2[["d_stage1"]] <- NULL

  if(latent.utility.targeted){
    # replace with working utility
    X.t_stage2[,"D"] <- rnorm(nrow(X.t_stage2))
    df.stage2$X[,"D"] <- X.t_stage2[,"D"]
  }

  df.stage1 <- df
  names(df.stage1)[names(df.stage1) == "X_stage1"] <- "X"
  names(df.stage1)[names(df.stage1) == "d_stage1"] <- "d"
  X.t_stage1 <- cbind(df.stage1$X, t = df.stage1$timeidx)
  df.stage1[["X_stage2"]] <- NULL
  df.stage1[["d_stage2"]] <- NULL

  mcmc_risk <- matrix(nrow = mcmc.opt$chain.length, ncol = length(df$y))
  Y <- matrix(nrow = length(df$y), ncol = mcmc.opt$chain.length)

  #progress bar
  if(progress.bar){
    pb <- utils::txtProgressBar(min = 0,
                                max = mcmc.opt$chain.length,
                                char = "=",
                                style = 3,
                                width = 30)
  }

  time <- system.time({

    for(i in 1:mcmc.opt$chain.length){

      # Step 1: Sampling of D.star (latent utility) ----------------------------

      if(i == 1){
         betat_stage1 <- matrix(rnorm(df$d_stage1*df$Tmax), nrow = df$Tmax, ncol = df$d_stage1)
      }
      b.t_stage1 <- cbind(c(t(betat_stage1)), rep(1:df$Tmax, each = df$d_stage1))
      colnames(b.t_stage1) <- c("b", "t")
      eta_stage1 <- lapply(1:df$Tmax, FUN = function(t){
        as.matrix(X.t_stage1[X.t_stage1[,"t"] == t, -ncol(X.t_stage1)]) %*%
          b.t_stage1[b.t_stage1[,"t"] == t, -ncol(b.t_stage1)]
      })
      eta_stage1 <- do.call("rbind", eta_stage1)

      if(i == 1){
        betat_stage2 <- matrix(rnorm(df$d_stage2*df$Tmax), nrow = df$Tmax, ncol = df$d_stage2)
      }
      reff.t_stage2 <- cbind(reff_stage2, t = X.t_stage2[,"t"])
      b.t_stage2 <- cbind(c(t(betat_stage2)), rep(1:df$Tmax, each = df$d_stage2))
      colnames(b.t_stage2) <- c("b", "t")
      eta_stage2 <- lapply(1:df$Tmax, FUN = function(t){
        as.matrix(X.t_stage2[X.t_stage2[,"t"] == t, -ncol(X.t_stage2)]) %*%
          b.t_stage2[b.t_stage2[,"t"] == t, -ncol(b.t_stage2)] +
          reff.t_stage2[reff.t_stage2[,"t"]==t,-ncol(reff.t_stage2)]
      })
      eta_stage2 <- do.call("rbind", eta_stage2)

      mu <- eta_stage1 + (rho/sqrt(sigma2v)) * (df$y - eta_stage2)
      D.star <- numeric(length(df$D))
      idx.zero <- df$D == 0
      idx.one <- df$D == 1
      D.star[idx.zero] <- truncnorm::rtruncnorm(n = sum(idx.zero), b = 0,
                                                mean = mu[idx.zero,], sd = sqrt(1-rho^2))
      D.star[idx.one] <- truncnorm::rtruncnorm(n = sum(idx.one), a = 0,
                                               mean = mu[idx.one,], sd = sqrt(1-rho^2))
      if(latent.utility.targeted){
        # replace with working utility
        X.t_stage2[,"D"] <- D.star
        df.stage2$X[,"D"] <- D.star
      }

      # Step 2: Sampling 2nd stage effects -------------------------------------

      y.tilde <- df$y - rho * sqrt(sigma2v) * (D.star - eta_stage1)
      s2.adj <- sigma2v * (1 - rho^2) # still homoscedastic

      # Step R

      ybeta_stage2 <- y.tilde - reff_stage2
      stepR.out_stage2 <- stepR(response = ybeta_stage2,
                               df = df.stage2,
                               prior.reg = prior.reg_stage2,
                               sigma2v = s2.adj,
                               alpha = alpha_stage2,
                               estimation = "Normal",
                               mcmc.opt = mcmc.opt,
                               i = i)
      betat_stage2 <- stepR.out_stage2$betat
      alpha_stage2 <- stepR.out_stage2$alpha
      prior.reg_stage2 <- stepR.out_stage2$prior.reg

      # Step F

      linpred_stage2 <- construct.lp(X = df.stage2$X,
                                     Time = df.stage2$Tmax,
                                     timeidx = df.stage2$timeidx,
                                     betat = betat_stage2)

      if(random.effects){

        res.y_stage2 <- y.tilde - linpred_stage2
        stepF.out_stage2 <- stepF(response = res.y_stage2,
                                  df = df.stage2,
                                  sigma2v = s2.adj,
                                  lambda = lambda_stage2,
                                  alpha_lambda = alpha_lambda_stage2,
                                  prior.load = prior.load_stage2,
                                  estimation = "Normal",
                                  mcmc.opt = mcmc.opt)
        fi_stage2 <- stepF.out_stage2$fi
        lambda_stage2 <- stepF.out_stage2$lambda
        alpha_lambda_stage2 <- stepF.out_stage2$alpha_lambda
        prior.load_stage2 <- stepF.out_stage2$prior.load
        fv_stage2 <- rep(fi_stage2, df.stage2$Tmax)
        if(i>mcmc.opt$burnin & i%%mcmc.opt$thin==0){
          f_mat_stage2[fi.count_stage2,] <- fi_stage2
          fi.count_stage2 <- fi.count_stage2+1
          f_sum_stage2 <- f_sum_stage2+fi_stage2
        }
        if(!tv.load_stage2){
          reff_stage2 <- lambda_stage2*fv_stage2
        } else{
          reff_stage2 <- c(t(matrix(lambda_stage2, ncol=df.stage2$n, nrow=df.stage2$Tmax)))*fv_stage2
        }

      }

      # Step 3: Sampling 1st stage effects -------------------------------------

      # Step R (there is NO Step F in 1st stage equation!)

      stepR.out_stage1 <- stepR(response = D.star,
                                df = df.stage1,
                                prior.reg = prior.reg_stage1,
                                sigma2v = 1,
                                alpha = alpha_stage1,
                                estimation = "Normal",
                                mcmc.opt = mcmc.opt,
                                i = i)
      betat_stage1 <- stepR.out_stage1$betat
      alpha_stage1 <- stepR.out_stage1$alpha
      prior.reg_stage1 <- stepR.out_stage1$prior.reg

      linpred_stage1 <- construct.lp(X = df.stage1$X,
                                     Time = df.stage1$Tmax,
                                     timeidx = df.stage1$timeidx,
                                     betat = betat_stage1)

      # Step 4: Joint sampling of sigma2 and rho based on 2D-Slice-Sampling ----

      cov.par <- slice_IV_2D(response = df$y - linpred_stage2 - reff_stage2,
                             residual.stage1 = D.star - linpred_stage1,
                             sigma2 = sigma2v,
                             rho = rho,
                             prior.var_stage2 = prior.var_stage2,
                             prior.rho = prior.rho)
      sigma2v <- cov.par[["sigma2"]]
      rho <- cov.par[["rho"]]

      # Step 5: Data augmentation in case of missing response data -------------

      df$y[miss] <- StepAugment(eta.miss = c(linpred_stage2)[miss] + reff_stage2[miss],
                                model = "Gaussian",
                                sigma2 = sigma2v)

      # Returning --------------------------------------------------------------

      # STAGE 2

      if(prior.reg_stage2$type %in% c("rw1", "rw2")){ # shrinkage

        if(prior.reg_stage2$TG && !prior.reg_stage2$TG.alternative){ # triple Gamma

          res.i_stage2 <- c(i,
                           betat_stage2,
                           alpha_stage2[1:df.stage2$d],
                           alpha_stage2[(df.stage2$d+1):(2*df.stage2$d)],
                           prior.reg_stage2$tau,
                           prior.reg_stage2$xi,
                           prior.reg_stage2$a.tau,
                           prior.reg_stage2$kappa.tau,
                           prior.reg_stage2$a.xi,
                           prior.reg_stage2$kappa.xi,
                           prior.reg_stage2$c.tau,
                           prior.reg_stage2$kappa.tau.check,
                           prior.reg_stage2$c.xi,
                           prior.reg_stage2$kappa.xi.check,
                           sigma2v,
                           lambda_stage2)

        } else if(prior.reg_stage2$TG && prior.reg_stage2$TG.alternative){

          res.i_stage2 <- c(i,
                           betat_stage2,
                           alpha_stage2[1:df.stage2$d],
                           alpha_stage2[(df.stage2$d+1):(2*df.stage2$d)],
                           prior.reg_stage2$tau,
                           prior.reg_stage2$xi,
                           prior.reg_stage2$a.tau,
                           prior.reg_stage2$a.xi,
                           prior.reg_stage2$c.tau,
                           prior.reg_stage2$c.xi,
                           prior.reg_stage2$chi.tau.j,
                           prior.reg_stage2$chi.xi.j,
                           sigma2v,
                           lambda_stage2)

        } else{ # double Gamma

          res.i_stage2 <- c(i,
                           betat_stage2,
                           alpha_stage2[1:df.stage2$d],
                           alpha_stage2[(df.stage2$d+1):(2*df.stage2$d)],
                           prior.reg_stage2$tau,
                           prior.reg_stage2$xi,
                           prior.reg_stage2$a.tau,
                           prior.reg_stage2$kappa.tau,
                           prior.reg_stage2$a.xi,
                           prior.reg_stage2$kappa.xi,
                           sigma2v,
                           lambda_stage2
          )
        }

      } else{ # independence prior

        res.i_stage2 <- c(i,
                         betat_stage2,
                         sigma2v,
                         lambda_stage2
        )

      }

      if(prior.load_stage2$type %in% c("rw1", "rw2")){

        res.i_stage2 <- c(res.i_stage2,
                         m_lambda = alpha_lambda_stage2[1],
                         psi = alpha_lambda_stage2[2],
                         phi = prior.load_stage2$phi,
                         zeta = prior.load_stage2$zeta,
                         a.phi = prior.load_stage2$a.phi,
                         kappa.phi = prior.load_stage2$kappa.phi,
                         a.zeta = prior.load_stage2$a.zeta,
                         kappa.zeta = prior.load_stage2$kappa.zeta)

      }

      res_frame_stage2[i,] <- res.i_stage2

      # STAGE 1

      if(prior.reg_stage1$type %in% c("rw1", "rw2")){ # shrinkage

        if(prior.reg_stage1$TG && !prior.reg_stage1$TG.alternative){ # triple Gamma

          res.i_stage1 <- c(i,
                        betat_stage1,
                        alpha_stage1[1:df.stage1$d],
                        alpha_stage1[(df.stage1$d+1):(2*df.stage1$d)],
                        prior.reg_stage1$tau,
                        prior.reg_stage1$xi,
                        prior.reg_stage1$a.tau,
                        prior.reg_stage1$kappa.tau,
                        prior.reg_stage1$a.xi,
                        prior.reg_stage1$kappa.xi,
                        prior.reg_stage1$c.tau,
                        prior.reg_stage1$kappa.tau.check,
                        prior.reg_stage1$c.xi,
                        prior.reg_stage1$kappa.xi.check)

        } else if(prior.reg_stage1$TG && prior.reg_stage1$TG.alternative){

          res.i_stage1 <- c(i,
                        betat_stage1,
                        alpha_stage1[1:df.stage1$d],
                        alpha_stage1[(df.stage1$d+1):(2*df.stage1$d)],
                        prior.reg_stage1$tau,
                        prior.reg_stage1$xi,
                        prior.reg_stage1$a.tau,
                        prior.reg_stage1$a.xi,
                        prior.reg_stage1$c.tau,
                        prior.reg_stage1$c.xi,
                        prior.reg_stage1$chi.tau.j,
                        prior.reg_stage1$chi.xi.j)

        } else{ # double Gamma

          res.i_stage1 <- c(i,
                        betat_stage1,
                        alpha_stage1[1:df.stage1$d],
                        alpha_stage1[(df.stage1$d+1):(2*df.stage1$d)],
                        prior.reg_stage1$tau,
                        prior.reg_stage1$xi,
                        prior.reg_stage1$a.tau,
                        prior.reg_stage1$kappa.tau,
                        prior.reg_stage1$a.xi,
                        prior.reg_stage1$kappa.xi
          )
        }

      } else{ # independence prior

        res.i_stage1 <- c(i,
                      betat_stage1
        )

      }

      res_frame_stage1[i,] <- res.i_stage1

      rho_vec[i] <- rho
      D_save[i,] <- D.star

      Y[,i] <- df$y # important for missings and computation of WAIC
      if(progress.bar) utils::setTxtProgressBar(pb, i) # tracking progress


    } # for-loop

  }) # system.time

  if(progress.bar) close(pb)
  #print time
  cat(paste("MCMC sampling finished in", round(time[3]), "seconds."))

  # Setting Up Return Object ---------------------------------------------------

  nmc <- (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin

  # second stage model (Normal)

  res_stage2 <- res_frame_stage2[res_frame_stage2[,"SimNr"] > mcmc.opt$burnin,]
  res_stage2 <- res_stage2[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin),]
  res_mcmc_stage2 <- coda::mcmc(data = res_stage2[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  res_stage2[,startsWith(colnames(res_stage2), "theta")] <- abs(res_stage2[,startsWith(colnames(res_stage2), "theta")])
  colnames(res_stage2)[startsWith(colnames(res_stage2), "theta")] <- paste0("abs(theta", 1:df.stage2$d, ")")
  res_stage2[,startsWith(colnames(res_stage2), "psi")] <- abs(res_stage2[,startsWith(colnames(res_stage2), "psi")])
  colnames(res_stage2)[startsWith(colnames(res_stage2), "psi")] <- "abs(psi)"
  res_stage2 <- coda::mcmc(data = res_stage2[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  hpint_stage2 <- coda::HPDinterval(res_stage2, prob = HPD.coverage)
  mcmcsummary_stage2 <- cbind(hpint_stage2[,"lower"],
                             c(colMeans(res_stage2)),
                             c(apply(res_stage2,2,median)),
                             hpint_stage2[,"upper"],
                             c(apply(res_stage2,2,sd))
  )
  colnames(mcmcsummary_stage2) <- c("LO", "mean","median","UP","sd")
  fmean_stage2 <- f_sum_stage2/nmc

  # first stage model (Probit)

  res_stage1 <- res_frame_stage1[res_frame_stage1[,"SimNr"] > mcmc.opt$burnin,]
  res_stage1 <- res_stage1[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin),]
  res_mcmc_stage1 <- coda::mcmc(data = res_stage1[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  res_stage1[,startsWith(colnames(res_stage1), "theta")] <- abs(res_stage1[,startsWith(colnames(res_stage1), "theta")])
  colnames(res_stage1)[startsWith(colnames(res_stage1), "theta")] <- paste0("abs(theta", 1:df.stage1$d, ")")
  res_stage1 <- coda::mcmc(data = res_stage1[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  hpint_stage1 <- coda::HPDinterval(res_stage1, prob = HPD.coverage)
  mcmcsummary_stage1 <- cbind(hpint_stage1[,"lower"],
                          c(colMeans(res_stage1)),
                          c(apply(res_stage1,2,median)),
                          hpint_stage1[,"upper"],
                          c(apply(res_stage1,2,sd))
  )
  colnames(mcmcsummary_stage1) <- c("LO", "mean","median","UP","sd")

  # correlation of errors (degree of endogeneity)

  rho_vec <- rho_vec[(mcmc.opt$burnin+1):mcmc.opt$chain.length]
  rho_vec <- coda::mcmc(data = rho_vec[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin)])
  rho_mcmc <- coda::mcmc(data = rho_vec, start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  rho_interval <- coda::HPDinterval(rho_vec, prob = HPD.coverage)
  mcmcsummary_rho <- matrix(c(rho_interval[1], mean(rho_mcmc), median(rho_mcmc), rho_interval[2], sd(rho_mcmc)), nrow = 1)
  colnames(mcmcsummary_rho) <- c("LO", "mean","median","UP","sd")
  rownames(mcmcsummary_rho) <- "rho"

  # latent utility D.star
  D_save <- D_save[(mcmc.opt$burnin+1):mcmc.opt$chain.length,]
  D_save <- coda::mcmc(data = D_save[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin),])
  D_mcmc <- coda::mcmc(data = D_save, start = mcmc.opt$burnin+1, end = mcmc.opt$chain.length, thin = mcmc.opt$thin)

  # computing acceptance rates of Metropolis-based parameters
  acceptance.rates_stage1 <- matrix(nrow = 1, ncol = 4)
  if(prior.reg_stage1$type != "ind"){
    acceptance.rates_stage1[,1] <- accept.rate(accept = prior.reg_stage1$a.xi.accept, mcmc.opt = mcmc.opt)
    acceptance.rates_stage1[,2] <- accept.rate(accept = prior.reg_stage1$a.tau.accept, mcmc.opt = mcmc.opt)
    if(prior.reg_stage1$TG){
      acceptance.rates_stage1[,3] <- accept.rate(accept = prior.reg_stage1$c.xi.accept, mcmc.opt = mcmc.opt)
      acceptance.rates_stage1[,4] <- accept.rate(accept = prior.reg_stage1$c.tau.accept, mcmc.opt = mcmc.opt)
    }
  }
  acceptance.rates_stage2 <- matrix(nrow = 1, ncol = 4)
  if(prior.reg_stage2$type != "ind"){
    acceptance.rates_stage2 <- matrix(nrow = 1, ncol = 4)
    acceptance.rates_stage2[,1] <- accept.rate(accept = prior.reg_stage2$a.xi.accept, mcmc.opt = mcmc.opt)
    acceptance.rates_stage2[,2] <- accept.rate(accept = prior.reg_stage2$a.tau.accept, mcmc.opt = mcmc.opt)
    if(prior.reg_stage2$TG){
      acceptance.rates_stage2[,3] <- accept.rate(accept = prior.reg_stage2$c.xi.accept, mcmc.opt = mcmc.opt)
      acceptance.rates_stage2[,4] <- accept.rate(accept = prior.reg_stage2$c.tau.accept, mcmc.opt = mcmc.opt)
    }
  }
  acceptance.rates <- cbind(acceptance.rates_stage1, acceptance.rates_stage2)
  colnames(acceptance.rates) <- c("a.xi (stage1)", "a.tau (stage1)", "c.xi (stage1)", "c.tau (stage1)",
                                  "a.xi (stage2)", "a.tau (stage2)", "c.xi (stage2)", "c.tau (stage2)")
  if(sum(is.na(acceptance.rates)) == 8) acceptance.rates <- NULL

  # return
  df$y[miss] <- NA
  ret <- list(data = df,
              Y = Y,
              mcmc_stage1 = res_mcmc_stage1, mcmc_stage2 = res_mcmc_stage2, mcmc_rho = rho_mcmc,
              posterior_stage1 = mcmcsummary_stage1, posterior_stage2 = mcmcsummary_stage2,
              posterior_rho = mcmcsummary_rho,
              fmcmc_stage2 = f_mat_stage2[,1:df$n], fmean_stage2 = fmean_stage2,
              D_mcmc = D_mcmc,
              model = "Probit-IV",
              acceptance.rates = acceptance.rates,
              HPD.coverage = HPD.coverage,
              runtime = paste("Total Runtime for Bayesian Probit Instrumental Variable Model:",
                              round(time[3], 3), "seconds"),
              Treatment.Variable = Treatment.Variable)
  if(sum(miss) == 0) ret$Y <- NULL

  return(ret)

}


