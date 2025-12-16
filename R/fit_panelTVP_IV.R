fit_panelTVP_IV <- function(formula_stage1,
                            formula_stage2,
                            data = data,
                            prior.reg_stage1,
                            prior.reg_stage2,
                            prior.load_stage2,
                            prior.var_stage2,
                            prior.rho,
                            mcmc.opt,
                            HPD.coverage,
                            random.effects,
                            progress.bar){

  # Initialization

  if(prior.load_stage2$type == "cps"){
    tv.load_stage2 = FALSE
  } else{
    tv.load_stage2 = TRUE
  }

  resp <- all.vars(formula_stage2)[1]
  miss <- ifelse(is.na(data[,resp]), TRUE, FALSE)
  N.miss <- sum(miss)
  data[miss, resp] <- rnorm(n = N.miss)

  if(any(colnames(data) == "t")) dat <- data[,names(data) != "t"]
  if(any(colnames(dat) == "id")) dat <- dat[,names(dat) != "id"]

  mf_stage1 <- model.frame(formula = formula_stage1, data = dat, drop.unused.levels = TRUE)
  mt_stage1 <- attr(mf_stage1, "terms")
  x_stage1 <- model.matrix(mt_stage1, mf_stage1)
  mf_stage2 <- model.frame(formula = formula_stage2, data = dat, drop.unused.levels = TRUE)
  mt_stage2 <- attr(mf_stage2, "terms")
  x_stage2 <- model.matrix(mt_stage2, mf_stage2)

  D <- model.response(mf_stage1)

  tind <- data$t
  Tmax <- max(tind)
  id <- data$id
  y <- model.response(mf_stage2)
  y <- as.numeric(y) # to be safe
  y <- y[order(tind, id)]
  x_stage1 <- x_stage1[order(tind, id), , drop = FALSE]
  x_stage2 <- x_stage2[order(tind, id), , drop = FALSE]
  df <- data.frame(tind, id)
  df <- df[order(tind, id),]

  df <- list(y = y, D = D, X_stage1 = x_stage1, X_stage2 = x_stage2,
             Tmax = Tmax, n = length(y)/max(tind),
             size = length(y), d_stage1 = ncol(x_stage1), d_stage2 = ncol(x_stage2),
             timeidx = df$tind, idx = df$id)

  # sigma2 only for second stage!
  sigma2 <- 1
  sigma2v <- sigma2

  # rho (correlation of errors between stages)
  rho <- 0

  alpha_stage1 <- rnorm(df$d_stage1*2)
  alpha_stage2 <- rnorm(df$d_stage2*2)

  fi_stage2 <- rep(0, df$n)
  fv_stage2 <- rep(fi_stage2, df$Tmax)

  if(!tv.load_stage2){
    lambda_stage2 <- 0
    reff_stage2 <- lambda_stage2*fv_stage2
  } else{
    lambda_stage2 <- rep(0, df$Tmax)
    reff_stage2 <- c(t(matrix(lambda_stage2, ncol=df$n, nrow=df$Tmax)))*fv_stage2
  }

  if(random.effects){
    alpha_lambda_stage2 <- matrix(c(1.2,0.5))
  } else{
    alpha_lambda_stage2 <- matrix(c(0,0))
  }

  if(tv.load_stage2 & length(prior.load_stage2$L0) == 1){
    prior.load_stage2$L0 <- rep(prior.load_stage2$L0, Tmax)
  }

  prior.reg_stage1$tau <- rep(10, df$d_stage1)
  prior.reg_stage1$xi <- rep(10, df$d_stage1)
  prior.reg_stage2$tau <- rep(10, df$d_stage2)
  prior.reg_stage2$xi <- rep(10, df$d_stage2)
  prior.load_stage2$phi <- 1
  prior.load_stage2$zeta <- 1

  if(prior.reg_stage1$TG && !prior.reg_stage1$TG.alternative){
    prior.reg_stage1$kappa.tau.check <- rep(1, df$d_stage1)
    prior.reg_stage1$kappa.xi.check <- rep(1, df$d_stage1)
    prior.reg_stage1$ph.tau <- (2*prior.reg_stage1$c.tau)/(prior.reg_stage1$kappa.tau*prior.reg_stage1$a.tau)
    prior.reg_stage1$ph.xi <- (2*prior.reg_stage1$c.xi)/(prior.reg_stage1$kappa.xi*prior.reg_stage1$a.xi)
    prior.reg_stage1$tau.check <- (prior.reg_stage1$tau*prior.reg_stage1$kappa.tau.check)/prior.reg_stage1$ph.tau
    prior.reg_stage1$xi.check <- (prior.reg_stage1$xi*prior.reg_stage1$kappa.xi.check)/prior.reg_stage1$ph.xi
  }

  if(prior.reg_stage2$TG && !prior.reg_stage2$TG.alternative){
    prior.reg_stage2$kappa.tau.check <- rep(1, df$d_stage2)
    prior.reg_stage2$kappa.xi.check <- rep(1, df$d_stage2)
    prior.reg_stage2$ph.tau <- (2*prior.reg_stage2$c.tau)/(prior.reg_stage2$kappa.tau*prior.reg_stage2$a.tau)
    prior.reg_stage2$ph.xi <- (2*prior.reg_stage2$c.xi)/(prior.reg_stage2$kappa.xi*prior.reg_stage2$a.xi)
    prior.reg_stage2$tau.check <- (prior.reg_stage2$tau*prior.reg_stage2$kappa.tau.check)/prior.reg_stage2$ph.tau
    prior.reg_stage2$xi.check <- (prior.reg_stage2$xi*prior.reg_stage2$kappa.xi.check)/prior.reg_stage2$ph.xi
  }

  #-----------------------------------------------------------------------------
  # create return matrix for the MCMC samples

  if(prior.reg_stage1$type == "ind"){
    namesbetat_stage1 <- unlist(lapply(1:df$d_stage1, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat_stage1 <- unlist(lapply(1:df$d_stage1, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta_stage1 <- paste0("beta", 1:df$d_stage1)
    namestheta_stage1 <- paste0("theta", 1:df$d_stage1)
    namestau_stage1 <- paste0("tau2", 1:df$d_stage1)
    namesxi_stage1 <- paste0("xi2", 1:df$d_stage1)
  }
  if(prior.reg_stage2$type == "ind"){
    namesbetat_stage2 <- unlist(lapply(1:df$d_stage2, function(x) paste0("beta_t",x,1:df$Tmax)))
  } else{
    namesbetat_stage2 <- unlist(lapply(1:df$d_stage2, function(x) paste0("beta_t",x,1:df$Tmax)))
    namesbeta_stage2 <- paste0("beta", 1:df$d_stage2)
    namestheta_stage2 <- paste0("theta", 1:df$d_stage2)
    namestau_stage2 <- paste0("tau2", 1:df$d_stage2)
    namesxi_stage2 <- paste0("xi2", 1:df$d_stage2)
  }

  if(prior.reg_stage1$type == "ind"){
    cnames_stage1 <- c("SimNr", namesbetat_stage1)
  } else if(prior.reg_stage1$TG && !prior.reg_stage1$TG.alternative){ # original Triple Gamma
    kappa.tau.j_stage1 <- paste0("kappa.tau_", 1:df$d_stage1)
    kappa.xi.j_stage1 <- paste0("kappa.xi_", 1:df$d_stage1)
    cnames_stage1 <- c("SimNr",namesbetat_stage1, namesbeta_stage1, namestheta_stage1,
                   namestau_stage1, namesxi_stage1, c("a.tau","kappa.tau","a.xi","kappa.xi",
                                              "c.tau", kappa.tau.j_stage1, "c.xi", kappa.xi.j_stage1))
  } else if(prior.reg_stage1$TG && prior.reg_stage1$TG.alternative){ # alternative Triple Gamma
    chi.tau.j_stage1 <- paste0("chi.tau_", 1:df$d_stage1)
    chi.xi.j_stage1 <- paste0("chi.xi_", 1:df$d_stage1)
    cnames_stage1 <- c("SimNr", namesbetat_stage1, namesbeta_stage1, namestheta_stage1,
                   namestau_stage1, namesxi_stage1, c("a.tau", "a.xi", "c.tau", "c.xi",
                                              chi.tau.j_stage1, chi.xi.j_stage1))
  } else{ # double Gamma
    cnames_stage1 <- c("SimNr",namesbetat_stage1, namesbeta_stage1, namestheta_stage1,
                   namestau_stage1, namesxi_stage1, c("a.tau","kappa.tau","a.xi","kappa.xi"))
  }
  col_res_stage1 <- length(cnames_stage1)
  res_frame_stage1 <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res_stage1)
  colnames(res_frame_stage1) <- cnames_stage1

  namessgma2_stage2 <- "sigma2" # only on stage 2!
  if(!tv.load_stage2){
    nameslambdat_stage2 <-"lambda_t"
  } else{
    nameslambdat_stage2 <- paste0("lambda_t",1:df$Tmax)
  }
  if(prior.reg_stage2$type == "ind"){
    cnames_stage2 <- c("SimNr", namesbetat_stage2, namessgma2_stage2, nameslambdat_stage2)
  } else if(prior.reg_stage2$TG && !prior.reg_stage2$TG.alternative){ # original Triple Gamma
    kappa.tau.j_stage2 <- paste0("kappa.tau_", 1:df$d_stage2)
    kappa.xi.j_stage2 <- paste0("kappa.xi_", 1:df$d_stage2)
    cnames_stage2 <- c("SimNr",namesbetat_stage2, namesbeta_stage2, namestheta_stage2,
                      namestau_stage2, namesxi_stage2, c("a.tau","kappa.tau","a.xi","kappa.xi",
                                                       "c.tau", kappa.tau.j_stage2, "c.xi", kappa.xi.j_stage2),
                      namessgma2_stage2, nameslambdat_stage2)
  } else if(prior.reg_stage2$TG && prior.reg_stage2$TG.alternative){ # alternative Triple Gamma
    chi.tau.j_stage2 <- paste0("chi.tau_", 1:df$d_stage2)
    chi.xi.j_stage2 <- paste0("chi.xi_", 1:df$d_stage2)
    cnames_stage2 <- c("SimNr", namesbetat_stage2, namesbeta_stage2, namestheta_stage2,
                      namestau_stage2, namesxi_stage2, c("a.tau", "a.xi", "c.tau", "c.xi",
                                                       chi.tau.j_stage2, chi.xi.j_stage2),
                      namessgma2_stage2, nameslambdat_stage2)
  } else{ # double Gamma
    cnames_stage2 <- c("SimNr",namesbetat_stage2, namesbeta_stage2, namestheta_stage2,
                      namestau_stage2, namesxi_stage2, c("a.tau","kappa.tau","a.xi","kappa.xi"),
                      namessgma2_stage2, nameslambdat_stage2)
  }
  col_res_stage2 <- length(cnames_stage2)
  if(prior.load_stage2$type=="rw1" || prior.load_stage2$type=="rw2"){
    cnames_stage2 <- c(cnames_stage2,"lambda","psi","phi2","zeta2", "a.phi", "kappa.phi",
                      "a.zeta", "kappa.zeta")
    col_res_stage2 <- length(cnames_stage2)
  }
  res_frame_stage2 <- matrix(0, nrow = mcmc.opt$chain.length, ncol = col_res_stage2)
  colnames(res_frame_stage2) <- cnames_stage2

  f_sum_stage2 <- rep(0, df$n)
  f_mat_stage2 <- matrix(NA, nrow = (mcmc.opt$chain.length-mcmc.opt$burnin)/mcmc.opt$thin,
                        ncol = length(y))

  ## regression part
  if(prior.reg_stage1$type != "ind"){
    prior.reg_stage1$a.xi.accept <- c()
    prior.reg_stage1$a.xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
    prior.reg_stage1$a.tau.accept <- c()
    prior.reg_stage1$a.tau.accept[1] <- 1 # we let metropolis start in 2nd iteration
    if(prior.reg_stage1$TG){
      prior.reg_stage1$c.xi.accept <- c()
      prior.reg_stage1$c.xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
      prior.reg_stage1$c.tau.accept <- c()
      prior.reg_stage1$c.tau.accept[1] <- 1 # we let metropolis start in 2nd iteration
    }
  }
  if(prior.reg_stage2$type != "ind"){
    prior.reg_stage2$a.xi.accept <- c()
    prior.reg_stage2$a.xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
    prior.reg_stage2$a.tau.accept <- c()
    prior.reg_stage2$a.tau.accept[1] <- 1 # we let metropolis start in 2nd iteration
    if(prior.reg_stage2$TG){
      prior.reg_stage2$c.xi.accept <- c()
      prior.reg_stage2$c.xi.accept[1] <- 1 # we let metropolis start in 2nd iteration
      prior.reg_stage2$c.tau.accept <- c()
      prior.reg_stage2$c.tau.accept[1] <- 1 # we let metropolis start in 2nd iteration
    }
  }

  # Fitting the IV model

  Treatment.Variable <- as.character(formula_stage1[[2]])
  result <- InstrumentalTVP(df = df,
                            prior.reg_stage1 = prior.reg_stage1,
                            prior.reg_stage2 = prior.reg_stage2,
                            prior.var_stage2 = prior.var_stage2,
                            prior.load_stage2 = prior.load_stage2,
                            prior.rho = prior.rho,
                            mcmc.opt = mcmc.opt,
                            alpha_stage1 = alpha_stage1,
                            alpha_stage2 = alpha_stage2,
                            lambda_stage2 = lambda_stage2,
                            alpha_lambda_stage2 = alpha_lambda_stage2,
                            reff_stage2 = reff_stage2,
                            tv.load_stage2 = tv.load_stage2,
                            res_frame_stage1 = res_frame_stage1,
                            res_frame_stage2 = res_frame_stage2,
                            f_sum_stage2 = f_sum_stage2,
                            f_mat_stage2 = f_mat_stage2,
                            miss = miss,
                            sigma2v = sigma2v,
                            rho = rho,
                            HPD.coverage = HPD.coverage,
                            random.effects = random.effects,
                            progress.bar = progress.bar,
                            Treatment.Variable = Treatment.Variable)
  class(result) <- "panelTVP.IV"

  return(result)

}
