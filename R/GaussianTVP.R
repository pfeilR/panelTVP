GaussianTVP <- function(df,
                        prior.reg,
                        prior.var,
                        prior.load,
                        mcmc.opt,
                        sigma2v,
                        alpha,
                        lambda,
                        alpha_lambda,
                        reff,
                        C0,
                        tv.load,
                        res_frame,
                        f_sum,
                        f_mat,
                        HPD.coverage){

  #progress bar
  pb <- utils::txtProgressBar(min = 0,
                       max = mcmc.opt$chain.length,
                       char = "=",
                       style = 3,
                       width = 30)

  fi.count <- 1

  time <- system.time({

    for(i in 1:mcmc.opt$chain.length){

      # Step R

      ybeta <- df$y - reff
      stepR.out <- stepR(response = ybeta,
                         df = df,
                         prior.reg = prior.reg,
                         sigma2v = sigma2v,
                         alpha = alpha,
                         estimation = "Normal",
                         mcmc.opt = mcmc.opt,
                         i = i)
      betat <- stepR.out$betat
      alpha <- stepR.out$alpha
      prior.reg <- stepR.out$prior.reg

      # Step V

      linpred <- construct.lp(X = df$X,
                              Time = df$Tmax,
                              timeidx = df$timeidx,
                              betat = betat)
      res.y <- df$y - linpred
      stepV.out <- stepV(response = res.y-reff,
                         df = df,
                         prior.var = prior.var,
                         C0 = C0)
      sigma2v <- stepV.out$sigma2
      C0 <- stepV.out$C0
      prior.var <- stepV.out$prior.var

      # Step F

      stepF.out <- stepF(response = res.y,
                         df = df,
                         sigma2v = sigma2v,
                         lambda = lambda,
                         alpha_lambda = alpha_lambda,
                         prior.load = prior.load,
                         estimation = "Normal")
      fi <- stepF.out$fi
      lambda <- stepF.out$lambda
      alpha_lambda <- stepF.out$alpha_lambda
      prior.load <- stepF.out$prior.load
      fv <- rep(fi, df$Tmax)
      if(i>mcmc.opt$burnin & i%%mcmc.opt$thin==0){
        f_mat[fi.count,] <- fi
        fi.count <- fi.count+1
        f_sum <- f_sum+fi
      }
      if(!tv.load){
        reff <- lambda*fv
      }else{
        reff <- c(t(matrix(lambda, ncol=df$n, nrow=df$Tmax)))*fv
      }

      # Returning

      if(prior.reg$type %in% c("rw1", "rw2")){ # shrinkage

        res.i <- c(i,
                   betat,
                   alpha[1:df$d],
                   alpha[(df$d+1):(2*df$d)],
                   prior.reg$tau,
                   prior.reg$xi,
                   prior.reg$a_tau,
                   prior.reg$kappa_tau,
                   prior.reg$a_xi,
                   prior.reg$kappa_xi,
                   sigma2v,
                   lambda
        )

      } else{ # independence prior

        res.i <- c(i,
                   betat,
                   sigma2v,
                   lambda
        )

      }

      if(prior.load$type %in% c("rw1", "rw2")){

        res.i <- c(res.i,
                   m_lambda = alpha_lambda[1],
                   psi = alpha_lambda[2],
                   phi_lambda = prior.load$phi,
                   zeta_lambda = prior.load$zeta,
                   a_phi = prior.load$a_phi,
                   kappa_phi = prior.load$kappa_phi,
                   a_zeta = prior.load$a_zeta,
                   kappa_zeta = prior.load$kappa_zeta)

      }

      res_frame[i,] <- res.i

      utils::setTxtProgressBar(pb, i) # tracking progress

    } # for-loop

  }) # system.time

  close(pb)
  #print time
  print(paste("Algorithm took", time[3], "seconds"))

  #remove burnin
  res <- res_frame[res_frame[,"SimNr"] > mcmc.opt$burnin,]
  res <- res[seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin),]
  idx <- res[,"SimNr"]
  res_mcmc <- coda::mcmc(data = res[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)

  #apply the absolute value to theta for the summary
  res[,startsWith(colnames(res), "theta")] <- abs(res[,startsWith(colnames(res), "theta")])
  colnames(res)[startsWith(colnames(res), "theta")] <- paste0("abs(theta", 1:df$d, ")")

  # added on 05.08.24: apply the absolute value to psi for the summary
  res[,startsWith(colnames(res), "psi")] <- abs(res[,startsWith(colnames(res), "psi")])
  colnames(res)[startsWith(colnames(res), "psi")] <- "abs(psi)"

  # in summary report lambda with positive sign for lambda 1
  if(prior.load$type != "cps"){
    sign.lambda1 <- sign(res[, "lambda_t1"])
    ind <- startsWith(colnames(res), "lambda")
    res[,ind] <- sign.lambda1*res[,ind]
  } else{
    res[, "lambda_t"] <- abs(res[, "lambda_t"])
  }
  res <- coda::mcmc(data = res[,-1], start = mcmc.opt$burnin+1, thin = mcmc.opt$thin)
  hpint <- coda::HPDinterval(res, prob = HPD.coverage)

  #----------------------------------------------
  #create own summary statistics
  mcmcsummary <- cbind(hpint[,"lower"],
                       c(colMeans(res)),
                       c(apply(res,2,median)),
                       hpint[,"upper"],
                       c(apply(res,2,sd))
  )
  colnames(mcmcsummary) <- c("LO", "mean","median","UP","sd")

  #create return object:
  nmc <- (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin
  fmean <- f_sum/nmc

  # computing acceptance rates of Metropolis-based parameters
  acceptance.rates <- matrix(nrow = 1, ncol = 2)
  acceptance.rates[,1] <- accept.rate(accept = prior.reg$xi.accept, mcmc.opt = mcmc.opt)
  acceptance.rates[,2] <- accept.rate(accept = prior.reg$tau.accept, mcmc.opt = mcmc.opt)
  colnames(acceptance.rates) <- c("a_xi", "a_tau")

  # return
  ret <- list(data = df, mcmc = res_mcmc, posterior = mcmcsummary,
              fmean = fmean, fmcmc = f_mat[,1:df$n], model = "Gaussian", acceptance.rates = acceptance.rates,
              HPD.coverage = HPD.coverage,
              runtime = paste("Total Runtime for Bayesian Normal Model:", round(time[3], 3), "seconds"))

  return(ret)

} # end
