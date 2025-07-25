% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3_print.R
\name{print.panelTVP.Logit}
\alias{print.panelTVP.Logit}
\title{Print basic model information output for a \code{panelTVP.Logit} object}
\usage{
\method{print}{panelTVP.Logit}(x, ...)
}
\arguments{
\item{x}{an object of class \code{panelTVP.Logit}}

\item{...}{optional arguments passed to the function (those are ignored)}
}
\description{
This basic \code{print} method gives a general overview
of the information that is contained in an object of class \code{panelTVP.Logit}.
}
\examples{
# Printing object of class panelTVP.Logit
# NB: To reduce computational effort, we have drastically reduced the length
# of the Markov Chain. You should use a much longer chain in your applications.
sim.logit <- sim_panelTVP(n = 100,
                          Tmax = 4,
                          beta = c(1,0.5,0,0),
                          theta = c(0.8,0.5,0,0),
                          lambda = 1,
                          psi = 0.2,
                          model = "Logit")
res.logit <- panelTVP(y ~ W1 + W2 + W3,
                      data = sim.logit$observed,
                      mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
                      model = "Logit")
print(res.logit)
}
\author{
Roman Pfeiler, Helga Wagner
}
