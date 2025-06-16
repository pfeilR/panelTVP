stepV <- function(response,
                  df,
                  prior.var,
                  C0){

  sigma2 <- sample_Sigma(ycen = response, n = df$n, Time = df$Tmax,
                         C0 = C0, hyper_par.c0 = prior.var$c0)

  C0 <- sample_C0(prior.var, x = sigma2)

  return(list(sigma2 = sigma2, C0 = C0, prior.var = prior.var))

}

# Additional Functions called by stepV -----------------------------------------

sample_Sigma <- function(ycen, n, Time, C0, hyper_par.c0){

  a <- hyper_par.c0+(n*Time/2)
  b <- C0 + 0.5*crossprod(ycen, ycen)
  sigma2 <- 1/rgamma(1,a,b)

  return(sigma2)

}
sample_C0 <- function(hyper_par, x){

  a <- hyper_par$learn.C0.hyp$g0 + hyper_par$c0
  b <- hyper_par$learn.C0.hyp$G0 + (1/x)
  rgamma(1,a,b)

}
