fit_panelTVP_ZINB <- function(formula,
                              data,
                              prior.reg_nb,
                              prior.load_nb,
                              prior.reg_logit,
                              prior.load_logit,
                              mcmc.opt,
                              settings.NegBin,
                              HPD.coverage,
                              progress.bar){

  # Initialization

  if(prior.load_nb$type == "cps"){
    tv.load_nb = FALSE
  } else{
    tv.load_nb = TRUE
  }

  if(prior.load_logit$type == "cps"){
    tv.load_logit = FALSE
  } else{
    tv.load_logit = TRUE
  }

  resp <- all.vars(formula)[1]
  miss <- ifelse(is.na(data[,resp]), TRUE, FALSE)
  N.miss <- sum(miss)
  data$y[miss] <- MASS::rnegbin(N.miss, mu = 1, theta = 2)

  rhs <- formula[[3]]
  rhs_str <- paste(deparse(rhs), collapse = "")
  rhs.parts <- strsplit(rhs_str, "\\|")[[1]]
  rhs.parts <- trimws(rhs.parts)
  formula_nb <- as.formula(paste(deparse(formula[[2]]), "~", rhs.parts[1]))
  formula_logit  <- as.formula(paste(deparse(formula[[2]]), "~", rhs.parts[2]))
  mf_nb <- model.frame(formula = formula_nb, data = data, drop.unused.levels = TRUE)
  mt_nb <- attr(mf_nb, "terms")
  x_nb <- model.matrix(mt_nb, mf_nb)
  mf_logit <- model.frame(formula = formula_logit, data = data, drop.unused.levels = TRUE)
  mt_logit <- attr(mf_logit, "terms")
  x_logit <- model.matrix(mt_logit, mf_logit)
  y <- model.response(mf_nb)
  tind <- data$t
  Tmax <- max(tind)
  id <- data$id
  y <- y[order(tind, id)]
  x_nb <- x_nb[order(tind, id), , drop = FALSE]
  x_logit <- x_logit[order(tind, id), , drop = FALSE]
  df <- data.frame(tind, id)
  df <- df[order(tind, id),]
  df <- list(y = y, X_nb = x_nb, X_logit = x_logit,
             Tmax = Tmax, n = length(y)/max(tind), size = length(y),
             d_nb = ncol(x_nb), d_logit = ncol(x_logit),
             timeidx = df$tind, idx = df$id)

  alpha_nb <- rnorm(df$d_nb*2)
  alpha_logit <- rnorm(df$d_logit*2)

  fi_nb <- rep(0, df$n)
  fv_nb <- rep(fi_nb, df$Tmax)
  fi_logit <- rep(0, df$n)
  fv_logit <- rep(fi_logit, df$Tmax)

  if(!tv.load_nb){
    lambda_nb <- 0
    reff_nb <- lambda_nb*fv_nb
  } else{
    lambda_nb <- rep(0, df$Tmax)
    reff_nb <- c(t(matrix(lambda_nb, ncol=df$n, nrow=df$Tmax)))*fv_nb
  }
  if(!tv.load_logit){
    lambda_logit <- 0
    reff_logit <- lambda_logit*fv_logit
  } else{
    lambda_logit <- rep(0, df$Tmax)
    reff_logit <- c(t(matrix(lambda_logit, ncol=df$n, nrow=df$Tmax)))*fv_logit
  }

  alpha_lambda_nb <- matrix(c(1.2,0.5))
  alpha_lambda_logit <- matrix(c(1.2,0.5))

  if(tv.load_nb & length(prior.load_nb$L0) == 1){
    prior.load_nb$L0 <- rep(prior.load_nb$L0, Tmax)
  }
  if(tv.load_logit & length(prior.load_logit$L0) == 1){
    prior.load_logit$L0 <- rep(prior.load_logit$L0, Tmax)
  }

  prior.reg_nb$tau <- rep(10, df$d_nb)
  prior.reg_nb$xi <- rep(10, df$d_nb)
  prior.load_nb$phi <- 1
  prior.load_nb$zeta <- 1

  prior.reg_logit$tau <- rep(10, df$d_logit)
  prior.reg_logit$xi <- rep(10, df$d_logit)
  prior.load_logit$phi <- 1
  prior.load_logit$zeta <- 1

  #-----------------------------------------------------------------------------
  # create return matrix for the MCMC samples

  if(prior.reg_nb$type == "ind"){
    namesbetat_nb <- unlist(lapply(1:df$d_nb, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat_nb <- unlist(lapply(1:df$d_nb, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta_nb <- paste0("beta", 1:df$d_nb)
    namestheta_nb <- paste0("theta", 1:df$d_nb)
    namestau_nb <- paste0("tau2", 1:df$d_nb)
    namesxi_nb <- paste0("xi2", 1:df$d_nb)
  }
  if(prior.reg_logit$type == "ind"){
    namesbetat_logit <- unlist(lapply(1:df$d_logit, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat_logit <- unlist(lapply(1:df$d_logit, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta_logit <- paste0("beta", 1:df$d_logit)
    namestheta_logit <- paste0("theta", 1:df$d_logit)
    namestau_logit <- paste0("tau2", 1:df$d_logit)
    namesxi_logit <- paste0("xi2", 1:df$d_logit)
  }

  if(!tv.load_nb){
    nameslambdat_nb <-"lambda_t"
  } else{
    nameslambdat_nb <- paste0("lambda_t",1:df$Tmax)
  }
  if(prior.reg_nb$type == "ind"){
    cnames_nb <- c("SimNr", namesbetat_nb, nameslambdat_nb)
  } else{
    cnames_nb <- c("SimNr",namesbetat_nb, namesbeta_nb, namestheta_nb,
                   namestau_nb, namesxi_nb, c("a.tau","kappa.tau","a.xi","kappa.xi"),
                   nameslambdat_nb)
  }
  if(prior.load_nb$type=="rw1" | prior.load_nb$type=="rw2"){
    cnames_nb <- c(cnames_nb,"lambda","psi","phi2","zeta2", "a.phi", "kappa.phi",
                   "a.zeta", "kappa.zeta")
  }
  cnames_nb <- append(cnames_nb, "r")
  col_res_nb <- length(cnames_nb)
  res_frame_nb <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res_nb)
  colnames(res_frame_nb) <- cnames_nb
  if(!tv.load_logit){
    nameslambdat_logit <-"lambda_t"
  } else{
    nameslambdat_logit <- paste0("lambda_t",1:df$Tmax)
  }
  if(prior.reg_logit$type == "ind"){
    cnames_logit <- c("SimNr", namesbetat_logit, nameslambdat_logit)
  } else{
    cnames_logit <- c("SimNr",namesbetat_logit, namesbeta_logit, namestheta_logit,
                      namestau_logit, namesxi_logit, c("a.tau","kappa.tau","a.xi","kappa.xi"),
                      nameslambdat_logit)
  }
  col_res_logit <- length(cnames_logit)
  if(prior.load_logit$type=="rw1" | prior.load_logit$type=="rw2"){
    cnames_logit <- c(cnames_logit,"lambda","psi","phi2","zeta2", "a.phi", "kappa.phi",
                      "a.zeta", "kappa.zeta")
    col_res_logit <- length(cnames_logit)
  }
  res_frame_logit <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res_logit)
  colnames(res_frame_logit) <- cnames_logit

  f_sum_nb <- rep(0, df$n)
  f_mat_nb <- matrix(NA, nrow = (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin,
                     ncol = length(y))
  f_sum_logit <- rep(0, df$n)
  f_mat_logit <- matrix(NA, nrow = (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin,
                        ncol = length(y))

  # initialize acceptance rate for dispersion parameter r
  settings.NegBin$r.accept <- c()

  ## regression part
  prior.reg_nb$xi.accept <- c()
  prior.reg_nb$xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg_nb$tau.accept <- c()
  prior.reg_nb$tau.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg_logit$xi.accept <- c()
  prior.reg_logit$xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
  prior.reg_logit$tau.accept <- c()
  prior.reg_logit$tau.accept[1] <- 1 # we let metropolis start in 2nd iteration

  # Fitting the ZINB model

  result <- zinbTVP(df = df,
                    prior.reg_nb = prior.reg_nb,
                    prior.reg_logit = prior.reg_logit,
                    prior.load_nb = prior.load_nb,
                    prior.load_logit = prior.load_logit,
                    mcmc.opt = mcmc.opt,
                    settings.NegBin = settings.NegBin,
                    alpha_nb = alpha_nb,
                    alpha_logit = alpha_logit,
                    lambda_nb = lambda_nb,
                    lambda_logit = lambda_logit,
                    alpha_lambda_nb = alpha_lambda_nb,
                    alpha_lambda_logit = alpha_lambda_logit,
                    reff_nb = reff_nb,
                    reff_logit = reff_logit,
                    tv.load_nb = tv.load_nb,
                    tv.load_logit = tv.load_logit,
                    res_frame_nb = res_frame_nb,
                    res_frame_logit = res_frame_logit,
                    f_sum_nb = f_sum_nb,
                    f_sum_logit = f_sum_logit,
                    f_mat_nb = f_mat_nb,
                    f_mat_logit = f_mat_logit,
                    miss = miss,
                    HPD.coverage = HPD.coverage,
                    progress.bar = progress.bar)
  class(result) <- "panelTVP.ZINB"

  return(result)

}

