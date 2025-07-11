
<!-- README.md is generated from README.Rmd. Please edit that file -->

# panelTVP

Hi folks and welcome to the package panelTVP!

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
#> Warning: Paket 'LaplacesDemon' wurde unter R Version 4.4.3 erstellt
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
#>             y         W1          W2 t id
#> 1   2.1821176 -0.5098477  1.76766171 1  1
#> 2   2.6398232  1.0831044  2.25618252 1  2
#> 3   7.2587612  0.8584773 -0.77560624 1  3
#> 4   4.3149185  1.6750408  0.25700549 1  4
#> 5  -0.1749802  0.4164253  0.21220292 1  5
#> 6   1.6287371  0.5274689  0.10420287 1  6
#> 7   3.5242557 -1.4124874 -0.05310962 1  7
#> 8   3.7894223 -0.8682197  1.04394507 1  8
#> 9   3.3283360  1.4774228  1.13492082 1  9
#> 10  4.9456265 -0.5645181 -1.08109683 1 10
```

The time varying regression effects and factor loadings can be accessed
as follows:

``` r
d$beta
#>          [,1]      [,2] [,3]
#> [1,] 3.045087 0.4649952    0
#> [2,] 3.807460 0.5885984    0
#> [3,] 3.500623 0.8339338    0
#> [4,] 4.179776 0.2488185    0
d$lambda
#>           [,1]
#> [1,] 0.9079866
#> [2,] 1.0911503
#> [3,] 1.1232735
#> [4,] 0.8847742
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
                model = "Gaussian",
                mcmc.opt = list(chain.length = 500, burnin = 100, thin = 2, asis = TRUE))
#>   |                                      |                              |   0%  |                                      |                              |   1%  |                                      |                              |   2%  |                                      |=                             |   2%  |                                      |=                             |   3%  |                                      |=                             |   4%  |                                      |=                             |   5%  |                                      |==                            |   5%  |                                      |==                            |   6%  |                                      |==                            |   7%  |                                      |==                            |   8%  |                                      |===                           |   8%  |                                      |===                           |   9%  |                                      |===                           |  10%  |                                      |===                           |  11%  |                                      |===                           |  12%  |                                      |====                          |  12%  |                                      |====                          |  13%  |                                      |====                          |  14%  |                                      |====                          |  15%  |                                      |=====                         |  15%  |                                      |=====                         |  16%  |                                      |=====                         |  17%  |                                      |=====                         |  18%  |                                      |======                        |  18%  |                                      |======                        |  19%  |                                      |======                        |  20%  |                                      |======                        |  21%  |                                      |======                        |  22%  |                                      |=======                       |  22%  |                                      |=======                       |  23%  |                                      |=======                       |  24%  |                                      |=======                       |  25%  |                                      |========                      |  25%  |                                      |========                      |  26%  |                                      |========                      |  27%  |                                      |========                      |  28%  |                                      |=========                     |  28%  |                                      |=========                     |  29%  |                                      |=========                     |  30%  |                                      |=========                     |  31%  |                                      |=========                     |  32%  |                                      |==========                    |  32%  |                                      |==========                    |  33%  |                                      |==========                    |  34%  |                                      |==========                    |  35%  |                                      |===========                   |  35%  |                                      |===========                   |  36%  |                                      |===========                   |  37%  |                                      |===========                   |  38%  |                                      |============                  |  38%  |                                      |============                  |  39%  |                                      |============                  |  40%  |                                      |============                  |  41%  |                                      |============                  |  42%  |                                      |=============                 |  42%  |                                      |=============                 |  43%  |                                      |=============                 |  44%  |                                      |=============                 |  45%  |                                      |==============                |  45%  |                                      |==============                |  46%  |                                      |==============                |  47%  |                                      |==============                |  48%  |                                      |===============               |  48%  |                                      |===============               |  49%  |                                      |===============               |  50%  |                                      |===============               |  51%  |                                      |===============               |  52%  |                                      |================              |  52%  |                                      |================              |  53%  |                                      |================              |  54%  |                                      |================              |  55%  |                                      |=================             |  55%  |                                      |=================             |  56%  |                                      |=================             |  57%  |                                      |=================             |  58%  |                                      |==================            |  58%  |                                      |==================            |  59%  |                                      |==================            |  60%  |                                      |==================            |  61%  |                                      |==================            |  62%  |                                      |===================           |  62%  |                                      |===================           |  63%  |                                      |===================           |  64%  |                                      |===================           |  65%  |                                      |====================          |  65%  |                                      |====================          |  66%  |                                      |====================          |  67%  |                                      |====================          |  68%  |                                      |=====================         |  68%  |                                      |=====================         |  69%  |                                      |=====================         |  70%  |                                      |=====================         |  71%  |                                      |=====================         |  72%  |                                      |======================        |  72%  |                                      |======================        |  73%  |                                      |======================        |  74%  |                                      |======================        |  75%  |                                      |=======================       |  75%  |                                      |=======================       |  76%  |                                      |=======================       |  77%  |                                      |=======================       |  78%  |                                      |========================      |  78%  |                                      |========================      |  79%  |                                      |========================      |  80%  |                                      |========================      |  81%  |                                      |========================      |  82%  |                                      |=========================     |  82%  |                                      |=========================     |  83%  |                                      |=========================     |  84%  |                                      |=========================     |  85%  |                                      |==========================    |  85%  |                                      |==========================    |  86%  |                                      |==========================    |  87%  |                                      |==========================    |  88%  |                                      |===========================   |  88%  |                                      |===========================   |  89%  |                                      |===========================   |  90%  |                                      |===========================   |  91%  |                                      |===========================   |  92%  |                                      |============================  |  92%  |                                      |============================  |  93%  |                                      |============================  |  94%  |                                      |============================  |  95%  |                                      |============================= |  95%  |                                      |============================= |  96%  |                                      |============================= |  97%  |                                      |============================= |  98%  |                                      |==============================|  98%  |                                      |==============================|  99%  |                                      |==============================| 100%
#> [1] "Algorithm took 3.57 seconds"
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
#> t=1      2.9104         2.9990      3.1003 0.0502
#> t=2      3.7887         3.9033      4.0092 0.0548
#> t=3      3.2724         3.3877      3.5107 0.0580
#> t=4      4.0798         4.1839      4.2822 0.0525
#> 
#> ==================================================
#>                  Estimates for W1
#> ==================================================
#> 
#>     Lower (HPD) Posterior Mean Upper (HPD)     SD
#> t=1      0.3136         0.4260      0.5108 0.0526
#> t=2      0.5345         0.6428      0.7836 0.0617
#> t=3      0.6623         0.7903      0.8835 0.0575
#> t=4      0.2112         0.3277      0.4232 0.0559
#> 
#> ==================================================
#>                  Estimates for W2
#> ==================================================
#> 
#>     Lower (HPD) Posterior Mean Upper (HPD)     SD
#> t=1     -0.0192         0.0145      0.0718 0.0242
#> t=2     -0.0267         0.0200      0.0810 0.0319
#> t=3     -0.0326         0.0143      0.0691 0.0248
#> t=4     -0.0205         0.0167      0.0771 0.0271
#> 
#> ==================================================
#>                  Estimates for λ
#> ==================================================
#> 
#>     Lower (HPD) Posterior Mean Upper (HPD)     SD
#> t=1      0.8101         0.9463      1.0680 0.0674
#> t=2      0.8889         1.0270      1.1449 0.0680
#> t=3      0.8852         0.9967      1.1108 0.0572
#> t=4      0.8183         0.9535      1.0645 0.0626
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
#>             Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)      2.9104         2.9990      3.1003 0.0502
#> W1               0.3136         0.4260      0.5108 0.0526
#> W2              -0.0192         0.0145      0.0718 0.0242
#> λ                0.8101         0.9463      1.0680 0.0674
#> 
#> ==================================================
#>             Estimates for Timepoint 2
#> ==================================================
#> 
#>             Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)      3.7887         3.9033      4.0092 0.0548
#> W1               0.5345         0.6428      0.7836 0.0617
#> W2              -0.0267         0.0200      0.0810 0.0319
#> λ                0.8889         1.0270      1.1449 0.0680
#> 
#> ==================================================
#>             Estimates for Timepoint 3
#> ==================================================
#> 
#>             Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)      3.2724         3.3877      3.5107 0.0580
#> W1               0.6623         0.7903      0.8835 0.0575
#> W2              -0.0326         0.0143      0.0691 0.0248
#> λ                0.8852         0.9967      1.1108 0.0572
#> 
#> ==================================================
#>             Estimates for Timepoint 4
#> ==================================================
#> 
#>             Lower (HPD) Posterior Mean Upper (HPD)     SD
#> (Intercept)      4.0798         4.1839      4.2822 0.0525
#> W1               0.2112         0.3277      0.4232 0.0559
#> W2              -0.0205         0.0167      0.0771 0.0271
#> λ                0.8183         0.9535      1.0645 0.0626
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
a creating a first overview of the package. For a more thorough
introduction to the package visit the packages’ Vignette under **INSERT
LINK TO VIGNETTE**.
