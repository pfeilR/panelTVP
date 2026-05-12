StepAugment <- function(eta.miss, model, sigma2 = NULL, r = NULL, risk.miss = NULL){

  m <- length(eta.miss)
  y.miss <- vector("numeric", m)

  if(model == "Gaussian"){
    y.miss <- rnorm(m, mean = eta.miss, sd = sqrt(sigma2))
  }

  if(model == "Probit"){
    y.miss <- rbinom(m, size = 1, prob = pnorm(eta.miss))
  }

  if(model == "Logit"){
    y.miss <- rbinom(m, size = 1, prob = plogis(eta.miss))
  }

  if(model == "NegBin"){
    y.miss <- MASS::rnegbin(m, mu = r * exp(eta.miss), theta = r)
  }

  if(model == "ZINB"){
    risi <- risk.miss == TRUE
    y.miss <- vector("numeric", m)
    y.miss[risi] <- MASS::rnegbin(sum(risi), mu = r * exp(eta.miss[risi]), theta = r)
  }

  return(y.miss)

}
