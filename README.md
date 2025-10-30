
<!-- README.md is generated from README.Rmd. Please edit that file -->

# panelTVP

Hello and welcome to the package panelTVP!

This package provides functions for fitting Bayesian regression models
for panel data with time-varying parameters. The functions in this
package were created by Roman Pfeiler and Helga Wagner. Maintainer of
the package is Roman Pfeiler (<roman.pfeiler@jku.at>).

## Installation

You can install the development version of panelTVP from Roman Pfeiler’s
personal GitHub page [GitHub](https://github.com/pfeilR/)

``` r
install_github("pfeilR/panelTVP")
```

## Introduction to the package

The main goal of this package is to fit flexible Bayesian regression
models for panel data in the time-varying parameter framework. The model
is flexible in the sense that for each covariate and the
subject-specific random effect, time-varying effects are estimated. For
the random effects, a factor model specification is used, where
time-varying factor loadings are multiplied with subject-specific factor
scores. This strategy allows estimation of time-varying random effects.
To avoid overfitting, hierarchical shrinkage priors are used to identify
whether an effect is time-varying, time-invariant or zero. The following
models are implemented:

- Gaussian model for normally distributed responses
- Logit / Probit model for binary responses
- Negative Binomial model for overdispersed count data
- Zero-Inflated Negative Binomial model for zero-inflated and
  overdispersed count data

## Small case-study: Gaussian response data

In the following, we show how to fit a Bayesian time-varying parameter
model for panel data, where the conditional distribution of the response
variable is assumed to be Gaussian. For simplicity, the statistical
analysis is based on simulated data.

Data can be simulated by using the function sim_panelTVP() from the
panelTVP package as follows:

``` r
devtools::load_all()
#> ℹ Loading panelTVP
d <- sim_panelTVP(n = 1000,
                  Tmax = 4,
                  model = "Gaussian",
                  beta = c(4,1,0),
                  theta = c(1,0.5,0),
                  lambda = 1,
                  psi = 0.2,
                  sigma2 = 3)
```

Object **d** contains the generated dataset that includes the covariates
and columns with time and subject indices, respectively:

``` r
head(d$observed, 10)
#>            y         W1         W2 t id
#> 1   5.196982 -1.4937710 -0.6106686 1  1
#> 2   7.858228 -1.2946407 -0.2107574 1  2
#> 3   5.324331  0.9601771  0.8036318 1  3
#> 4   9.733167  0.4050892 -1.8623255 1  4
#> 5   6.893841 -0.5909410 -2.2632236 1  5
#> 6  10.698213  0.2038481 -0.7180575 1  6
#> 7   7.412297  0.2297553 -1.3912966 1  7
#> 8   7.451582  0.5505888 -0.7580365 1  8
#> 9   8.172961  0.5287267 -1.3559472 1  9
#> 10  9.955344  2.7738593 -1.1137174 1 10
```

The time varying regression effects and factor loadings can be accessed
as follows:

``` r
d$beta
#>          [,1]      [,2] [,3]
#> [1,] 6.920883 0.6775745    0
#> [2,] 7.164630 0.1583311    0
#> [3,] 6.874853 0.1648405    0
#> [4,] 6.436847 0.2691158    0
d$lambda
#>           [,1]
#> [1,] 0.9423072
#> [2,] 0.7772188
#> [3,] 1.3045569
#> [4,] 0.9192251
```

It should be noted that the first covariate effect corresponds to the
global intercept.

For fitting a Bayesian time-varying parameter model to the panel data,
the main function panelTVP() is used. The function contains a lot of
arguments, which are explained in greater depth in the corresponding
help page. Moreover, only a very small number of draws in MCMC
estimation are used here as this is only an illustrative example. In
applications, you should opt for a longer Markov Chain.

``` r
res <- panelTVP(y ~ W1 + W2,
                data = d$observed,
                id = d$observed$id,
                t = d$observed$t,
                model = "Gaussian",
                mcmc.opt = list(chain.length = 500, burnin = 100, thin = 2, asis = TRUE))
#> MCMC sampling finished in 3 seconds. Computing WAIC ... Analysis finished!
```

A summary table of the estimated coefficients can easily be created by
calling the corresponding S3-function:

``` r
summary(res)
#> 
#> ------------------------------------------------------------------------------
#> Posterior Summary of the Bayesian Normal Model with Time-Varying Coefficients:
#> ------------------------------------------------------------------------------
#> ==================================================
#>             Estimates for (Intercept)
#> ==================================================
#> 
#>     Lower (HPD) Posterior Mean Upper (HPD)     SD
#> t=1      6.7955         6.8893      6.9892 0.0494
#> t=2      7.0720         7.1628      7.2740 0.0539
#> t=3      6.8330         6.9307      7.0266 0.0522
#> t=4      6.2139         6.3170      6.4249 0.0547
#> 
#> ==================================================
#>                  Estimates for W1
#> ==================================================
#> 
#>     Lower (HPD) Posterior Mean Upper (HPD)     SD
#> t=1      0.6554         0.7643      0.9241 0.0644
#> t=2      0.1922         0.2862      0.4103 0.0577
#> t=3      0.0662         0.1850      0.3015 0.0610
#> t=4      0.1613         0.2819      0.4092 0.0615
#> 
#> ==================================================
#>                  Estimates for W2
#> ==================================================
#> 
#>     Lower (HPD) Posterior Mean Upper (HPD)     SD
#> t=1     -0.0707         0.0089      0.0780 0.0367
#> t=2     -0.0727        -0.0014      0.0654 0.0357
#> t=3     -0.0837        -0.0068      0.0675 0.0383
#> t=4     -0.1385        -0.0311      0.0471 0.0484
#> 
#> ==================================================
#>            Estimates for Factor Loading
#> ==================================================
#> 
#>     Lower (HPD) Posterior Mean Upper (HPD)     SD
#> t=1      0.8309         0.9517      1.0966 0.0696
#> t=2      0.6741         0.7891      0.9821 0.0791
#> t=3      1.1657         1.3127      1.4480 0.0782
#> t=4      0.6928         0.8111      0.9552 0.0705
```

Alternatively, results can also be ordered by time point:

``` r
summary(res, by = "timepoint")
#> 
#> ------------------------------------------------------------------------------
#> Posterior Summary of the Bayesian Normal Model with Time-Varying Coefficients:
#> ------------------------------------------------------------------------------
#> ==================================================
#>             Estimates for Timepoint 1
#> ==================================================
#> 
#>                Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)         6.7955         6.8893      6.9892 0.0494
#> W1                  0.6554         0.7643      0.9241 0.0644
#> W2                 -0.0707         0.0089      0.0780 0.0367
#> Factor Loading      0.8309         0.9517      1.0966 0.0696
#> 
#> ==================================================
#>             Estimates for Timepoint 2
#> ==================================================
#> 
#>                Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)         7.0720         7.1628      7.2740 0.0539
#> W1                  0.1922         0.2862      0.4103 0.0577
#> W2                 -0.0727        -0.0014      0.0654 0.0357
#> Factor Loading      0.6741         0.7891      0.9821 0.0791
#> 
#> ==================================================
#>             Estimates for Timepoint 3
#> ==================================================
#> 
#>                Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)         6.8330         6.9307      7.0266 0.0522
#> W1                  0.0662         0.1850      0.3015 0.0610
#> W2                 -0.0837        -0.0068      0.0675 0.0383
#> Factor Loading      1.1657         1.3127      1.4480 0.0782
#> 
#> ==================================================
#>             Estimates for Timepoint 4
#> ==================================================
#> 
#>                Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)         6.2139         6.3170      6.4249 0.0547
#> W1                  0.1613         0.2819      0.4092 0.0615
#> W2                 -0.1385        -0.0311      0.0471 0.0484
#> Factor Loading      0.6928         0.8111      0.9552 0.0705
```

Covariate effect plots - that include the estimated factor loadings as
well - are created as follows:

``` r
plot(res, nplots = 2)
```

<img src="man/figures/README-unnamed-chunk-8-1.png" width="100%" />

    #> The plots are out there, hit [Return] to see ...

<img src="man/figures/README-unnamed-chunk-8-2.png" width="100%" />

As we can see, the insignificant effect of covariate **W2** is
effectively shrunken towards zero, while the other effects vary over
time - as we would expect based on how the data were generated.

## The end?

Note that this was just a very simple and quick example for the sake of
creating a first overview of the package. For a more thorough
introduction to the package visit the packages’ Vignette under **INSERT
LINK TO VIGNETTE**.
