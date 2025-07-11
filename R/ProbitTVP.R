# ProbitTVP - This function contains the full MCMC sampler for the Probit model.

# Last update: 06.04.25 (RP)

ProbitTVP <- function(df,
                      prior.reg,
                      prior.load,
                      mcmc.opt,
                      alpha,
                      lambda,
                      alpha_lambda,
                      reff,
                      tv.load,
                      res_frame,
                      f_sum,
                      f_mat,
                      miss,
                      HPD.coverage){

  X.t <- cbind(df$X, t = df$timeidx)
  fi.count <- 1
  Y <- matrix(nrow = length(df$y), ncol = mcmc.opt$chain.length)

  #progress bar
  pb <- utils::txtProgressBar(min = 0,
                       max = mcmc.opt$chain.length,
                       char = "=",
                       style = 3,
                       width = 30)

  time <- system.time({

    for(i in 1:mcmc.opt$chain.length){

      # Step U

      if(i == 1){
        betat <- matrix(rnorm(df$d*df$Tmax), nrow = df$Tmax, ncol = df$d)
      }
      reff.t <- cbind(reff, t = X.t[,"t"])
      b.t <- cbind(c(t(betat)), rep(1:df$Tmax, each = df$d))
      colnames(b.t) <- c("b", "t")
      eta <- lapply(1:df$Tmax, FUN = function(t){
        X.t[X.t[,"t"] == t, -ncol(X.t)] %*% b.t[b.t[,"t"] == t, -ncol(b.t)] +
          reff.t[reff.t[,"t"]==t,-ncol(reff.t)]
      })
      eta <- do.call("rbind", eta)
      z <- vector("numeric", length(df$y))
      y1 <- df$y == 1
      y0 <- df$y == 0
      z[y0] <- truncnorm::rtruncnorm(sum(y0), a = -Inf, b = 0, mean = eta[y0,], sd = 1)
      z[y1] <- truncnorm::rtruncnorm(sum(y1), a = 0, b = Inf, mean = eta[y1,], sd = 1)

      # Step R

      zbeta <- z - reff
      stepR.out <- stepR(response = zbeta,
                         df = df,
                         prior.reg = prior.reg,
                         sigma2v = 1,
                         alpha = alpha,
                         estimation = "Normal",
                         mcmc.opt = mcmc.opt,
                         i = i)
      # estimation = 'Normal' is correct for Probit as we set sigma2 = 1
      # such that we can then treat it as Normal with error variance 1
      betat <- stepR.out$betat
      alpha <- stepR.out$alpha
      prior.reg <- stepR.out$prior.reg

      # Step F

      linpred <- construct.lp(X = df$X,
                              Time = df$Tmax,
                              timeidx = df$timeidx,
                              betat = betat)
      res.z <- z - linpred
      stepF.out <- stepF(response = res.z,
                         df = df,
                         sigma2v = 1,
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

      # Step Augment (only in the presence of missings)
      df$y[miss] <- StepAugment(eta.miss = c(linpred)[miss] + reff[miss],
                                model = "Probit")

      # Returning

      if(prior.reg$type %in% c("rw1", "rw2")){ # shrinkage

        res.i <- c(i,
                   betat,
                   alpha[1:df$d],
                   alpha[(df$d+1):(2*df$d)],
                   prior.reg$tau,
                   prior.reg$xi,
                   prior.reg$a.tau,
                   prior.reg$kappa.tau,
                   prior.reg$a.xi,
                   prior.reg$kappa.xi,
                   0,
                   lambda
        )

      } else{ # independence prior

        res.i <- c(i,
                   betat,
                   0,
                   lambda
        )

      }

      if(prior.load$type %in% c("rw1", "rw2")){

        res.i <- c(res.i,
                   m_lambda = alpha_lambda[1],
                   psi = alpha_lambda[2],
                   phi = prior.load$phi,
                   zeta = prior.load$zeta,
                   a.phi = prior.load$a.phi,
                   kappa.phi = prior.load$kappa.phi,
                   a.zeta = prior.load$a.zeta,
                   kappa.zeta = prior.load$kappa.zeta)

      }

      res_frame[i,] <- res.i
      Y[,i] <- df$y
      utils::setTxtProgressBar(pb, i) # tracking progress

    } # for-loop

  }) # system.time

  close(pb)
  #print time
  print(paste("Algorithm took", time[3], "seconds"))

  # remove burnin
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
  Y <- Y[,seq(1, mcmc.opt$chain.length-mcmc.opt$burnin, by = mcmc.opt$thin)]

  # computing acceptance rates of Metropolis-based parameters
  acceptance.rates <- matrix(nrow = 1, ncol = 2)
  acceptance.rates[,1] <- accept.rate(accept = prior.reg$xi.accept, mcmc.opt = mcmc.opt)
  acceptance.rates[,2] <- accept.rate(accept = prior.reg$tau.accept, mcmc.opt = mcmc.opt)
  colnames(acceptance.rates) <- c("a.xi", "a.tau")

  # return
  df$y[miss] <- NA
  ret <- list(data = df, Y = Y, mcmc = res_mcmc[,colnames(res_mcmc) != "sigma2"],
              posterior = mcmcsummary[rownames(mcmcsummary) != "sigma2",],
              fmean = fmean,fmcmc = f_mat[,1:df$n], model = "Probit", acceptance.rates = acceptance.rates,
              HPD.coverage = HPD.coverage,
              runtime = paste("Total Runtime for Bayesian Probit Model:", round(time[3], 3), "seconds"))
  if(sum(miss) == 0) ret$Y <- NULL

  return(ret)

} # end
