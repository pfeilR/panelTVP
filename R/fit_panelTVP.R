fit_panelTVP <- function(formula,
                         data,
                         model,
                         prior.reg,
                         prior.var,
                         prior.load,
                         mcmc.opt,
                         settings.NegBin,
                         HPD.coverage,
                         random.effects,
                         progress.bar){

  # Initialization

  if(prior.load$type == "cps"){
    tv.load = FALSE
  } else{
    tv.load = TRUE
  }

  resp <- all.vars(formula)[1]
  miss <- ifelse(is.na(data[,resp]), TRUE, FALSE)
  N.miss <- sum(miss)
  if(model == "Gaussian") data$y[miss] <- rnorm(n = N.miss)
  if(model %in% c("Probit", "Logit")) data$y[miss] <- rbinom(n = N.miss, size = 1, prob = 0.5)
  if(model == "NegBin") data$y[miss] <- MASS::rnegbin(n = N.miss, mu = 1, theta = 1)

  mf <- model.frame(formula = formula, data = data, drop.unused.levels = TRUE)
  y <- model.response(mf)
  mt <- attr(mf, "terms")
  x <- model.matrix(mt, mf)

  tind <- data$t
  Tmax <- max(tind)
  id <- data$id

  y <- y[order(tind, id)]
  x <- x[order(tind, id), , drop = FALSE]
  df <- data.frame(tind, id)
  df <- df[order(tind, id),]

  df <- list(y = y, X = x, Tmax = Tmax,
             n = length(y)/max(tind),
             size = length(y), d = ncol(x),
             timeidx = df$tind, idx = df$id)

  alpha <- rnorm(df$d*2)
  sigma2 <- 1
  sigma2v <- sigma2

  C0 <- 1

  fi <- rep(0, df$n)
  fv <- rep(fi, df$Tmax)

  if(!tv.load){
    lambda <- 0
    reff <- lambda*fv
  } else{
    lambda <- rep(0, df$Tmax)
    reff <- c(t(matrix(lambda, ncol=df$n, nrow=df$Tmax)))*fv
  }

  alpha_lambda <- matrix(c(1.2,0.5))

  if(tv.load & length(prior.load$L0) == 1){
    prior.load$L0 <- rep(prior.load$L0, Tmax)
  }

  prior.reg$tau <- rep(10, df$d)
  prior.reg$xi <- rep(10, df$d)
  prior.load$phi <- 1
  prior.load$zeta <- 1

  # create return matrix for the MCMC samples

  if(prior.reg$type == "ind"){
    namesbetat <- unlist(lapply(1:df$d, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat <- unlist(lapply(1:df$d, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta <- paste0("beta", 1:df$d)
    namestheta <- paste0("theta", 1:df$d)
    namestau <- paste0("tau2", 1:df$d)
    namesxi <- paste0("xi2", 1:df$d)
  }
  namessgma2 <- "sigma2"
  if(!tv.load){
    nameslambdat <-"lambda_t"
  } else{
    nameslambdat <- paste0("lambda_t",1:df$Tmax)
  }
  if(prior.reg$type == "ind"){
    cnames <- c("SimNr", namesbetat, namessgma2, nameslambdat)
  } else{
    cnames <- c("SimNr",namesbetat, namesbeta, namestheta,
                namestau, namesxi, c("a.tau","kappa.tau","a.xi","kappa.xi"),
                namessgma2, nameslambdat)
  }
  if(prior.load$type=="rw1" | prior.load$type=="rw2"){
    cnames <- c(cnames,"lambda","psi","phi2","zeta2", "a.phi", "kappa.phi", "a.zeta", "kappa.zeta")
  }
  col_res <- length(cnames)
  res_frame <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res)
  colnames(res_frame) <- cnames

  f_sum <- rep(0, df$n)
  f_mat <- matrix(NA, nrow = (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin, ncol = length(y))

  # modification for Negative-Binomial model as we also want chain for r

  if(model == "NegBin"){
    cnames <- append(cnames, "r")
    col_res <- length(cnames)
    res_frame <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res)
    colnames(res_frame) <- cnames
  }

  # initialize acceptance rates

  ## r (negative binomial dispersion parameter)
  settings.NegBin$r.accept <- c()
  settings.NegBin$r.accept[1] <- 1 # we let metropolis start in 2nd iteration

  ## regression part
  prior.reg$xi.accept <- c()
  prior.reg$xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg$tau.accept <- c()
  prior.reg$tau.accept[1] <- 1 # we let metropolis start in 2nd iteration

  # fitting the model ----------------------------------------------------------

  if(model == "Gaussian"){
    result <- GaussianTVP(df = df,
                          prior.reg = prior.reg,
                          prior.var = prior.var,
                          prior.load = prior.load,
                          mcmc.opt = mcmc.opt,
                          sigma2v = sigma2v,
                          alpha = alpha,
                          lambda = lambda,
                          alpha_lambda = alpha_lambda,
                          reff = reff,
                          C0 = C0,
                          tv.load = tv.load,
                          res_frame = res_frame,
                          f_sum = f_sum,
                          f_mat = f_mat,
                          miss = miss,
                          HPD.coverage = HPD.coverage,
                          random.effects = random.effects,
                          progress.bar = progress.bar)
    class(result) <- "panelTVP.Gaussian"
  }

  if(model == "Probit"){
    result <- ProbitTVP(df = df,
                        prior.reg = prior.reg,
                        prior.load = prior.load,
                        mcmc.opt = mcmc.opt,
                        alpha = alpha,
                        lambda = lambda,
                        alpha_lambda = alpha_lambda,
                        reff = reff,
                        tv.load = tv.load,
                        res_frame = res_frame,
                        f_sum = f_sum,
                        f_mat = f_mat,
                        miss = miss,
                        HPD.coverage = HPD.coverage,
                        random.effects = random.effects,
                        progress.bar = progress.bar)
    class(result) <- "panelTVP.Probit"
  }

  if(model == "Logit"){
    result <- LogitTVP(df = df,
                       prior.reg = prior.reg,
                       prior.load = prior.load,
                       mcmc.opt = mcmc.opt,
                       alpha = alpha,
                       lambda = lambda,
                       alpha_lambda = alpha_lambda,
                       reff = reff,
                       tv.load = tv.load,
                       res_frame = res_frame,
                       f_sum = f_sum,
                       f_mat = f_mat,
                       miss = miss,
                       HPD.coverage = HPD.coverage,
                       random.effects = random.effects,
                       progress.bar = progress.bar)
    class(result) <- "panelTVP.Logit"
  }

  if(model == "NegBin"){
    result <- NegBinTVP(df = df,
                        prior.reg = prior.reg,
                        prior.load = prior.load,
                        mcmc.opt = mcmc.opt,
                        settings.NegBin = settings.NegBin,
                        alpha = alpha,
                        lambda = lambda,
                        alpha_lambda = alpha_lambda,
                        reff = reff,
                        tv.load = tv.load,
                        res_frame = res_frame,
                        f_sum = f_sum,
                        f_mat = f_mat,
                        miss = miss,
                        HPD.coverage = HPD.coverage,
                        random.effects = random.effects,
                        progress.bar = progress.bar)
    class(result) <- "panelTVP.NegBin"
  }

  return(result)

}


