zinbTVP <- function(df,
                    prior.reg_nb,
                    prior.reg_logit,
                    prior.load_nb,
                    prior.load_logit,
                    settings.NegBin,
                    mcmc.opt,
                    alpha_nb,
                    alpha_logit,
                    lambda_nb,
                    lambda_logit,
                    alpha_lambda_nb,
                    alpha_lambda_logit,
                    reff_nb,
                    reff_logit,
                    tv.load_nb,
                    tv.load_logit,
                    res_frame_nb,
                    res_frame_logit,
                    f_sum_nb,
                    f_sum_logit,
                    f_mat_nb,
                    f_mat_logit,
                    miss,
                    HPD.coverage,
                    random.effects,
                    progress.bar){

  r <- 1

  fi.count_nb <- 1
  fi.count_logit <- 1
  df.logit <- df
  names(df.logit)[names(df.logit) == "X_logit"] <- "X"
  names(df.logit)[names(df.logit) == "d_logit"] <- "d"
  X.t_logit <- cbind(df.logit$X, t = df.logit$timeidx)
  df.logit[["X_nb"]] <- NULL
  df.logit[["d_nb"]] <- NULL

  df.nb <- df
  names(df.nb)[names(df.nb) == "X_nb"] <- "X"
  names(df.nb)[names(df.nb) == "d_nb"] <- "d"
  X.t_nb <- cbind(df.nb$X, t = df.nb$timeidx)
  df.nb[["X_logit"]] <- NULL
  df.nb[["d_logit"]] <- NULL

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

      # Compute Linear Predictors for both Components --------------------------

      ## negative binomial

      if(i == 1){
        betat_nb <- matrix(rnorm(df$d_nb*df$Tmax), nrow = df$Tmax, ncol = df$d_nb)
      }
      reff.t_nb <- cbind(reff_nb, t = X.t_nb[,"t"])
      b.t_nb <- cbind(c(t(betat_nb)), rep(1:df$Tmax, each = df$d_nb))
      colnames(b.t_nb) <- c("b", "t")
      eta_nb <- lapply(1:df$Tmax, FUN = function(t){
        as.matrix(X.t_nb[X.t_nb[,"t"] == t, -ncol(X.t_nb)]) %*% b.t_nb[b.t_nb[,"t"] == t, -ncol(b.t_nb)] +
          reff.t_nb[reff.t_nb[,"t"]==t,-ncol(reff.t_nb)]
      })
      eta_nb <- do.call("rbind", eta_nb)

      ## logit

      if(i == 1){
        betat_logit <- matrix(rnorm(df$d_logit*df$Tmax), nrow = df$Tmax, ncol = df$d_logit)
      }
      reff.t_logit <- cbind(reff_logit, t = X.t_logit[,"t"])
      b.t_logit <- cbind(c(t(betat_logit)), rep(1:df$Tmax, each = df$d_logit))
      colnames(b.t_logit) <- c("b", "t")
      eta_logit <- lapply(1:df$Tmax, FUN = function(t){
        as.matrix(X.t_logit[X.t_logit[,"t"] == t, -ncol(X.t_logit)]) %*% b.t_logit[b.t_logit[,"t"] == t, -ncol(b.t_logit)] +
          reff.t_logit[reff.t_logit[,"t"]==t,-ncol(reff.t_logit)]
      })
      eta_logit <- do.call("rbind", eta_logit)

      # Sampling of latent at-risk indicators ----------------------------------

      risk <- stepRisk(y = df$y, miss = miss, eta_nb = eta_nb, eta_logit = eta_logit, r = r)
      mcmc_risk[i,] <- risk

      # Logit Component --------------------------------------------------------

      # Step U

      omega_logit <- pgdraw::pgdraw(b = 1, c = eta_logit)
      W.sparse_logit <- Matrix::Diagonal(n = length(omega_logit), x = omega_logit)
      W.dense_logit <- matrix(omega_logit, nrow = df$n, ncol = df$Tmax)
      z_logit <- (risk-1/2)/omega_logit

      # Step R

      zbeta_logit <- z_logit - reff_logit
      stepR.out_logit <- stepR(response = zbeta_logit,
                               df = df.logit,
                               prior.reg = prior.reg_logit,
                               W.sparse = W.sparse_logit,
                               alpha = alpha_logit,
                               estimation = "PG",
                               mcmc.opt = mcmc.opt,
                               i = i)
      betat_logit <- stepR.out_logit$betat
      alpha_logit <- stepR.out_logit$alpha
      prior.reg_logit <- stepR.out_logit$prior.reg

      # Step F

      linpred_logit <- construct.lp(X = df.logit$X,
                                    Time = df.logit$Tmax,
                                    timeidx = df.logit$timeidx,
                                    betat = betat_logit)

      if(random.effects){

        res.z_logit <- z_logit - linpred_logit
        stepF.out_logit <- stepF(response = res.z_logit,
                                 df = df.logit,
                                 W.sparse = W.sparse_logit,
                                 W.dense = W.dense_logit,
                                 lambda = lambda_logit,
                                 alpha_lambda = alpha_lambda_logit,
                                 prior.load = prior.load_logit,
                                 estimation = "PG")
        fi_logit <- stepF.out_logit$fi
        lambda_logit <- stepF.out_logit$lambda
        alpha_lambda_logit <- stepF.out_logit$alpha_lambda
        prior.load_logit <- stepF.out_logit$prior.load
        fv_logit <- rep(fi_logit, df.logit$Tmax)
        if(i>mcmc.opt$burnin & i%%mcmc.opt$thin==0){
          f_mat_logit[fi.count_logit,] <- fi_logit
          fi.count_logit <- fi.count_logit+1
          f_sum_logit <- f_sum_logit+fi_logit
        }
        if(!tv.load_logit){
          reff_logit <- lambda_logit*fv_logit
        } else{
          reff_logit <- c(t(matrix(lambda_logit, ncol=df.logit$n, nrow=df.logit$Tmax)))*fv_logit
        }

      }

      # Sampling of dispersion parameter r -------------------------------------

      y.risk <- df$y[risk]
      n.risk <- sum(risk)

      r.prev <- r
      sample.r.list <- NB.para(y = y.risk,
                               eta = matrix(eta_nb[risk,]),
                               r.old = r.prev,
                               sample.r = TRUE,
                               r.alpha = settings.NegBin$alpha.r,
                               r.beta = settings.NegBin$beta.r,
                               expansion.steps = settings.NegBin$expansion.steps,
                               width = settings.NegBin$width,
                               p.overrelax = settings.NegBin$p.overrelax,
                               accuracy.overrelax = settings.NegBin$accuracy.overrelax)
      r <- sample.r.list$r
      r <- r

      # Negative Binomial Component --------------------------------------------

      # Step U

      omega_nb <- vector(mode = "numeric", length = length(df$y))
      omega_nb[risk] <- efficient_PG_sampling(h = y.risk + r, z = c(eta_nb[risk,]))
      W.sparse_nb <- Matrix::Diagonal(n = length(omega_nb), x = omega_nb)
      W.dense_nb <- matrix(omega_nb, nrow = df$n, ncol = df$Tmax)
      z_nb <- vector(mode = "numeric", length = length(df$y))
      z_nb[risk] <- (df$y[risk]-r) / (2*omega_nb[risk])

      # Step R

      zbeta_nb <- z_nb - reff_nb
      stepR.out_nb <- stepR(response = zbeta_nb,
                            df = df.nb,
                            prior.reg = prior.reg_nb,
                            W.sparse = W.sparse_nb,
                            alpha = alpha_nb,
                            estimation = "PG",
                            mcmc.opt = mcmc.opt,
                            i = i)
      betat_nb <- stepR.out_nb$betat
      alpha_nb <- stepR.out_nb$alpha
      prior.reg_nb <- stepR.out_nb$prior.reg

      # Step F

      linpred_nb <- construct.lp(X = df.nb$X,
                                 Time = df.nb$Tmax,
                                 timeidx = df.nb$timeidx,
                                 betat = betat_nb)

      if(random.effects){

        res.z_nb <- z_nb - linpred_nb
        stepF.out_nb <- stepF(response = res.z_nb,
                              df = df.nb,
                              W.sparse = W.sparse_nb,
                              W.dense = W.dense_nb,
                              lambda = lambda_nb,
                              alpha_lambda = alpha_lambda_nb,
                              prior.load = prior.load_nb,
                              estimation = "PG")
        fi_nb <- stepF.out_nb$fi
        lambda_nb <- stepF.out_nb$lambda
        alpha_lambda_nb <- stepF.out_nb$alpha_lambda
        prior.load_nb <- stepF.out_nb$prior.load
        fv_nb <- rep(fi_nb, df.nb$Tmax)
        if(i>mcmc.opt$burnin & i%%mcmc.opt$thin==0){
          f_mat_nb[fi.count_nb,] <- fi_nb
          fi.count_nb <- fi.count_nb+1
          f_sum_nb <- f_sum_nb+fi_nb
        }
        if(!tv.load_nb){
          reff_nb <- lambda_nb*fv_nb
        }else{
          reff_nb <- c(t(matrix(lambda_nb, ncol=df.nb$n, nrow=df.nb$Tmax)))*fv_nb
        }

      }

      # Step Augment
      df$y[miss] <- StepAugment(eta.miss = c(linpred_nb)[miss] + reff_nb[miss],
                                model = "ZINB",
                                r = r,
                                risk.miss = risk[miss])

      # Returning --------------------------------------------------------------

      if(prior.reg_logit$type %in% c("rw1", "rw2")){ # shrinkage

        if(prior.reg_logit$TG && !prior.reg_logit$TG.alternative){ # triple Gamma

          res.i_logit <- c(i,
                           betat_logit,
                           alpha_logit[1:df.logit$d],
                           alpha_logit[(df.logit$d+1):(2*df.logit$d)],
                           prior.reg_logit$tau,
                           prior.reg_logit$xi,
                           prior.reg_logit$a.tau,
                           prior.reg_logit$kappa.tau,
                           prior.reg_logit$a.xi,
                           prior.reg_logit$kappa.xi,
                           prior.reg_logit$c.tau,
                           prior.reg_logit$kappa.tau.check,
                           prior.reg_logit$c.xi,
                           prior.reg_logit$kappa.xi.check,
                           lambda_logit)

        } else if(prior.reg_logit$TG && prior.reg_logit$TG.alternative){

          res.i_logit <- c(i,
                           betat_logit,
                           alpha_logit[1:df.logit$d],
                           alpha_logit[(df.logit$d+1):(2*df.logit$d)],
                           prior.reg_logit$tau,
                           prior.reg_logit$xi,
                           prior.reg_logit$a.tau,
                           prior.reg_logit$a.xi,
                           prior.reg_logit$c.tau,
                           prior.reg_logit$c.xi,
                           prior.reg_logit$chi.tau.j,
                           prior.reg_logit$chi.xi.j,
                           lambda_logit)

        } else{ # double Gamma

          res.i_logit <- c(i,
                           betat_logit,
                           alpha_logit[1:df.logit$d],
                           alpha_logit[(df.logit$d+1):(2*df.logit$d)],
                           prior.reg_logit$tau,
                           prior.reg_logit$xi,
                           prior.reg_logit$a.tau,
                           prior.reg_logit$kappa.tau,
                           prior.reg_logit$a.xi,
                           prior.reg_logit$kappa.xi,
                           lambda_logit
          )
        }

      } else{ # independence prior

        res.i_logit <- c(i,
                         betat_logit,
                         lambda_logit
        )

      }

      if(prior.load_logit$type %in% c("rw1", "rw2")){

        res.i_logit <- c(res.i_logit,
                         m_lambda = alpha_lambda_logit[1],
                         psi = alpha_lambda_logit[2],
                         phi = prior.load_logit$phi,
                         zeta = prior.load_logit$zeta,
                         a.phi = prior.load_logit$a.phi,
                         kappa.phi = prior.load_logit$kappa.phi,
                         a.zeta = prior.load_logit$a.zeta,
                         kappa.zeta = prior.load_logit$kappa.zeta)

      }

      res_frame_logit[i,] <- res.i_logit

      #

      if(prior.reg_nb$type %in% c("rw1", "rw2")){ # shrinkage

        if(prior.reg_nb$TG && !prior.reg_nb$TG.alternative){ # triple Gamma

          res.i_nb <- c(i,
                        betat_nb,
                        alpha_nb[1:df.nb$d],
                        alpha_nb[(df.nb$d+1):(2*df.nb$d)],
                        prior.reg_nb$tau,
                        prior.reg_nb$xi,
                        prior.reg_nb$a.tau,
                        prior.reg_nb$kappa.tau,
                        prior.reg_nb$a.xi,
                        prior.reg_nb$kappa.xi,
                        prior.reg_nb$c.tau,
                        prior.reg_nb$kappa.tau.check,
                        prior.reg_nb$c.xi,
                        prior.reg_nb$kappa.xi.check,
                        lambda_nb)

        } else if(prior.reg_nb$TG && prior.reg_nb$TG.alternative){

          res.i_nb <- c(i,
                        betat_nb,
                        alpha_nb[1:df.nb$d],
                        alpha_nb[(df.nb$d+1):(2*df.nb$d)],
                        prior.reg_nb$tau,
                        prior.reg_nb$xi,
                        prior.reg_nb$a.tau,
                        prior.reg_nb$a.xi,
                        prior.reg_nb$c.tau,
                        prior.reg_nb$c.xi,
                        prior.reg_nb$chi.tau.j,
                        prior.reg_nb$chi.xi.j,
                        lambda_nb)

        } else{ # double Gamma

          res.i_nb <- c(i,
                           betat_nb,
                           alpha_nb[1:df.nb$d],
                           alpha_nb[(df.nb$d+1):(2*df.nb$d)],
                           prior.reg_nb$tau,
                           prior.reg_nb$xi,
                           prior.reg_nb$a.tau,
                           prior.reg_nb$kappa.tau,
                           prior.reg_nb$a.xi,
                           prior.reg_nb$kappa.xi,
                           lambda_nb
          )
        }

      } else{ # independence prior

        res.i_nb <- c(i,
                      betat_nb,
                      lambda_nb
        )

      }

      if(prior.load_nb$type %in% c("rw1", "rw2")){

        res.i_nb <- c(res.i_nb,
                      m_lambda = alpha_lambda_nb[1],
                      psi = alpha_lambda_nb[2],
                      phi = prior.load_nb$phi,
                      zeta = prior.load_nb$zeta,
                      a.phi = prior.load_nb$a.phi,
                      kappa.phi = prior.load_nb$kappa.phi,
                      a.zeta = prior.load_nb$a.zeta,
                      kappa.zeta = prior.load_nb$kappa.zeta)

      }

      res_frame_nb[i,] <- c(res.i_nb, r)
      Y[,i] <- df$y # important for missings and computation of WAIC
      if(progress.bar) utils::setTxtProgressBar(pb, i) # tracking progress

    } # for-loop

  }) # system.time

  if(progress.bar) close(pb)
  #print time
  cat(paste("MCMC sampling finished in", round(time[3]), "seconds."))

  # Setting Up Return Object ---------------------------------------------------

  nmc <- (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin

  # logit component

  res_logit <- res_frame_logit[res_frame_logit[,"SimNr"] > mcmc.opt$burnin,]
  res_logit <- res_logit[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin),]
  res_mcmc_logit <- coda::mcmc(data = res_logit[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  res_logit[,startsWith(colnames(res_logit), "theta")] <- abs(res_logit[,startsWith(colnames(res_logit), "theta")])
  colnames(res_logit)[startsWith(colnames(res_logit), "theta")] <- paste0("abs(theta", 1:df.logit$d, ")")
  res_logit[,startsWith(colnames(res_logit), "psi")] <- abs(res_logit[,startsWith(colnames(res_logit), "psi")])
  colnames(res_logit)[startsWith(colnames(res_logit), "psi")] <- "abs(psi)"
  res_logit <- coda::mcmc(data = res_logit[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  hpint_logit <- coda::HPDinterval(res_logit, prob = HPD.coverage)
  mcmcsummary_logit <- cbind(hpint_logit[,"lower"],
                             c(colMeans(res_logit)),
                             c(apply(res_logit,2,median)),
                             hpint_logit[,"upper"],
                             c(apply(res_logit,2,sd))
  )
  colnames(mcmcsummary_logit) <- c("LO", "mean","median","UP","sd")
  fmean_logit <- f_sum_logit/nmc

  # negative binomial component

  res_nb <- res_frame_nb[res_frame_nb[,"SimNr"] > mcmc.opt$burnin,]
  res_nb <- res_nb[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin),]
  res_mcmc_nb <- coda::mcmc(data = res_nb[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  res_nb[,startsWith(colnames(res_nb), "theta")] <- abs(res_nb[,startsWith(colnames(res_nb), "theta")])
  colnames(res_nb)[startsWith(colnames(res_nb), "theta")] <- paste0("abs(theta", 1:df.nb$d, ")")
  res_nb[,startsWith(colnames(res_nb), "psi")] <- abs(res_nb[,startsWith(colnames(res_nb), "psi")])
  colnames(res_nb)[startsWith(colnames(res_nb), "psi")] <- "abs(psi)"
  res_nb <- coda::mcmc(data = res_nb[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  hpint_nb <- coda::HPDinterval(res_nb, prob = HPD.coverage)
  mcmcsummary_nb <- cbind(hpint_nb[,"lower"],
                          c(colMeans(res_nb)),
                          c(apply(res_nb,2,median)),
                          hpint_nb[,"upper"],
                          c(apply(res_nb,2,sd))
  )
  colnames(mcmcsummary_nb) <- c("LO", "mean","median","UP","sd")
  fmean_nb <- f_sum_nb/nmc

  # at-risk indicators
  mcmc_risk <- mcmc_risk[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin),]
  Y <- Y[,seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin)]

  # computing acceptance rates of Metropolis-based parameters
  acceptance.rates_logit <- matrix(nrow = 1, ncol = 4)
  if(prior.reg_logit$type != "ind" && !prior.reg_logit$TG.alternative){
    acceptance.rates_logit[,1] <- accept.rate(accept = prior.reg_logit$a.xi.accept, mcmc.opt = mcmc.opt)
    acceptance.rates_logit[,2] <- accept.rate(accept = prior.reg_logit$a.tau.accept, mcmc.opt = mcmc.opt)
    if(prior.reg_logit$TG){
      acceptance.rates_logit[,3] <- accept.rate(accept = prior.reg_logit$c.xi.accept, mcmc.opt = mcmc.opt)
      acceptance.rates_logit[,4] <- accept.rate(accept = prior.reg_logit$c.tau.accept, mcmc.opt = mcmc.opt)
    }
  }
  acceptance.rates_nb <- matrix(nrow = 1, ncol = 4)
  if(prior.reg_nb$type != "ind" && !prior.reg_nb$TG.alternative){
    acceptance.rates_nb <- matrix(nrow = 1, ncol = 4)
    acceptance.rates_nb[,1] <- accept.rate(accept = prior.reg_nb$a.xi.accept, mcmc.opt = mcmc.opt)
    acceptance.rates_nb[,2] <- accept.rate(accept = prior.reg_nb$a.tau.accept, mcmc.opt = mcmc.opt)
    if(prior.reg_nb$TG){
      acceptance.rates_nb[,3] <- accept.rate(accept = prior.reg_nb$c.xi.accept, mcmc.opt = mcmc.opt)
      acceptance.rates_nb[,4] <- accept.rate(accept = prior.reg_nb$c.tau.accept, mcmc.opt = mcmc.opt)
    }
  }
  acceptance.rates <- cbind(acceptance.rates_logit, acceptance.rates_nb)
  colnames(acceptance.rates) <- c("a.xi (logit)", "a.tau (logit)", "c.xi (logit)", "c.tau (logit)",
                                  "a.xi (nb)", "a.tau (nb)", "c.xi (nb)", "c.tau (nb)")
  if(sum(is.na(acceptance.rates)) == 8) acceptance.rates <- NULL

  # return
  df$y[miss] <- NA
  ret <- list(data = df,
              Y = Y,
              mcmc_logit = res_mcmc_logit, mcmc_nb = res_mcmc_nb,
              posterior_logit = mcmcsummary_logit, posterior_nb = mcmcsummary_nb,
              fmcmc_logit = f_mat_logit[,1:df$n], fmcmc_nb = f_mat_nb[,1:df$n],
              fmean_logit = fmean_logit, fmean_nb = fmean_nb,
              mcmc_risk = mcmc_risk,
              model = "ZINB",
              acceptance.rates = acceptance.rates,
              HPD.coverage = HPD.coverage,
              runtime = paste("Total Runtime for Bayesian Zero-Inflated Negative Binomial Model:",
                              round(time[3], 3), "seconds"))
  if(sum(miss) == 0) ret$Y <- NULL

  return(ret)

}
