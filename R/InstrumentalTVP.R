InstrumentalTVP <- function(df,
                            prior.reg_stage1,
                            prior.reg_stage2,
                            prior.load_stage1,
                            prior.load_stage2,
                            prior.var_stage2,
                            mcmc.opt,
                            alpha_stage1,
                            alpha_stage2,
                            lambda_stage1,
                            lambda_stage2,
                            alpha_lambda_stage1,
                            alpha_lambda_stage2,
                            reff_stage1,
                            reff_stage2,
                            tv.load_stage1,
                            tv.load_stage2,
                            res_frame_stage1,
                            res_frame_stage2,
                            f_sum_stage1,
                            f_sum_stage2,
                            f_mat_stage1,
                            f_mat_stage2,
                            miss,
                            sigma2v,
                            C0,
                            rho,
                            HPD.coverage,
                            random.effects_stage1,
                            random.effects_stage2,
                            progress.bar){

  fi.count_stage1 <- 1
  fi.count_stage2 <- 1

  df.stage2 <- df
  names(df.stage2)[names(df.stage2) == "X_stage2"] <- "X"
  names(df.stage2)[names(df.stage2) == "d_stage2"] <- "d"
  X.t_stage2 <- cbind(df.stage2$X, t = df.stage2$timeidx)
  df.stage2[["X_stage1"]] <- NULL
  df.stage2[["d_stage1"]] <- NULL

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
      reff.t_stage1 <- cbind(reff_stage1, t = X.t_stage1[,"t"])
      b.t_stage1 <- cbind(c(t(betat_stage1)), rep(1:df$Tmax, each = df$d_stage1))
      colnames(b.t_stage1) <- c("b", "t")
      eta_stage1 <- lapply(1:df$Tmax, FUN = function(t){
        as.matrix(X.t_stage1[X.t_stage1[,"t"] == t, -ncol(X.t_stage1)]) %*%
          b.t_stage1[b.t_stage1[,"t"] == t, -ncol(b.t_stage1)] +
          reff.t_stage1[reff.t_stage1[,"t"]==t,-ncol(reff.t_stage1)]
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
      D.star[idx.zero] <- truncnorm::rtruncnorm(n = sum(idx.zero), b = 0, mean = mu[idx.zero,], sd = sqrt(1-rho^2))
      D.star[idx.one] <- truncnorm::rtruncnorm(n = sum(idx.one), a = 0, mean = mu[idx.one,], sd = sqrt(1-rho^2))

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

      if(random.effects_stage2){

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

      # Step R

      ybeta_stage1 <- D.star - reff_stage1
      stepR.out_stage1 <- stepR(response = ybeta_stage1,
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

      # Step F

      linpred_stage1 <- construct.lp(X = df.stage1$X,
                                     Time = df.stage1$Tmax,
                                     timeidx = df.stage1$timeidx,
                                     betat = betat_stage1)

      if(random.effects_stage1){

        res.y_stage1 <- D.star - linpred_stage1
        stepF.out_stage1 <- stepF(response = res.y_stage1,
                                  df = df.stage1,
                                  sigma2v = 1,
                                  lambda = lambda_stage1,
                                  alpha_lambda = alpha_lambda_stage1,
                                  prior.load = prior.load_stage1,
                                  estimation = "Normal",
                                  mcmc.opt = mcmc.opt)
        fi_stage1 <- stepF.out_stage1$fi
        lambda_stage1 <- stepF.out_stage1$lambda
        alpha_lambda_stage1 <- stepF.out_stage1$alpha_lambda
        prior.load_stage1 <- stepF.out_stage1$prior.load
        fv_stage1 <- rep(fi_stage1, df.stage1$Tmax)
        if(i>mcmc.opt$burnin & i%%mcmc.opt$thin==0){
          f_mat_stage1[fi.count_stage1,] <- fi_stage1
          fi.count_stage1 <- fi.count_stage1+1
          f_sum_stage1 <- f_sum_stage1+fi_stage1
        }
        if(!tv.load_stage1){
          reff_stage1 <- lambda_stage1*fv_stage1
        } else{
          reff_stage1 <- c(t(matrix(lambda_stage1, ncol=df.stage1$n, nrow=df.stage1$Tmax)))*fv_stage1
        }

      }

      # Step 4: Sampling error variance ----------------------------------------

     # u <- df$y - linpred_stage2 - reff_stage2 - rho * (D.star - linpred_stage1 - reff_stage1)
    #  stepV.out <- stepV(response = u / sqrt(1-rho^2),
    #                     df = df,
     #                    prior.var = prior.var_stage2,
    #                     C0 = C0)
    #  sigma2v <- stepV.out$sigma2
    #  C0 <- stepV.out$C0
    #  prior.var_stage2 <- stepV.out$prior.var

      # Step 5: Sampling correlation between errors ----------------------------

      # Step 6: Data augmentation in case of missing response data -------------

      # Returning --------------------------------------------------------------

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

      #

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
                        prior.reg_stage1$kappa.xi.check,
                        lambda_stage1)

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
                        prior.reg_stage1$chi.xi.j,
                        lambda_stage1)

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
                        prior.reg_stage1$kappa.xi,
                        lambda_stage1
          )
        }

      } else{ # independence prior

        res.i_stage1 <- c(i,
                      betat_stage1,
                      lambda_stage1
        )

      }

      if(prior.load_stage1$type %in% c("rw1", "rw2")){

        res.i_stage1 <- c(res.i_stage1,
                      m_lambda = alpha_lambda_stage1[1],
                      psi = alpha_lambda_stage1[2],
                      phi = prior.load_stage1$phi,
                      zeta = prior.load_stage1$zeta,
                      a.phi = prior.load_stage1$a.phi,
                      kappa.phi = prior.load_stage1$kappa.phi,
                      a.zeta = prior.load_stage1$a.zeta,
                      kappa.zeta = prior.load_stage1$kappa.zeta)

      }

      res_frame_stage1[i,] <- res.i_stage1
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
  res_stage1[,startsWith(colnames(res_stage1), "psi")] <- abs(res_stage1[,startsWith(colnames(res_stage1), "psi")])
  colnames(res_stage1)[startsWith(colnames(res_stage1), "psi")] <- "abs(psi)"
  res_stage1 <- coda::mcmc(data = res_stage1[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  hpint_stage1 <- coda::HPDinterval(res_stage1, prob = HPD.coverage)
  mcmcsummary_stage1 <- cbind(hpint_stage1[,"lower"],
                          c(colMeans(res_stage1)),
                          c(apply(res_stage1,2,median)),
                          hpint_stage1[,"upper"],
                          c(apply(res_stage1,2,sd))
  )
  colnames(mcmcsummary_stage1) <- c("LO", "mean","median","UP","sd")
  fmean_stage1 <- f_sum_stage1/nmc

}


