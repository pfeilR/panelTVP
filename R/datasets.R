#' Marijuana dataset for binary panel models
#'
#' This dataset is a pre-processed subset of the National Longitudinal Survey of Youth 1997,
#' which is sponsored by the U.S. Bureau of Labor Statistics. This panel study provides a
#' rich data source on a wide range of socio-economic topics. This subsetted dataset serves
#' mainly as an application for binary time-varying parameter panel models. The goal is to
#' identify which covariates affect the probability of Marijuana consumption among
#' adolescents and young adults. The subjects are between 12 and 17 years old at study
#' begin, i.e., in the year 1997. 4609 subjects were measured at 7 consecutive time points.
#' Note that this is only a subset that was created for the demonstration of this very specific
#' drug application. The actual study comprises a total of 20 waves.
#'
#' \describe{
#'   \item{Bullied.Before.12}{was the respondent victim of bullying before age 12 (factor variable)}
#'   \item{Sex}{biological sex of the respondent (factor variable)}
#'   \item{Ethnicity}{ethnicity of the respondent (factor variable)}
#'   \item{Year}{year of survey round (numeric variable)}
#'   \item{Age}{age of the respondent (numeric variable)}
#'   \item{Residence}{does the respondent live in an urban or rural environment (factor variable)}
#'   \item{t}{time index (numeric variable)}
#'   \item{id}{respondent identifier (numeric variable); note that some id's
#'    are missing as some original respondents did not met the inclusion
#'    requirements}
#'   \item{Baseline.Age_c}{centered age at first panel wave, i.e.,
#'    a value of 0 means that the person was 12 years old at study begin (numeric variable)}
#'   \item{Used.Mari.Since.DLI}{has the respondent consumed Marijuana since the last interview /
#'    in the first wave the question was if the person has ever consumed Marijuana
#'    (numeric variable)}
#' }
#' @references Bureau of Labor Statistics, U.S. Department of Labor.
#'  National Longitudinal Survey of Youth 1997 cohort, 1997-2021 (rounds 1-20).
#'   Produced and distributed by the Center for Human Resource Research (CHRR),
#'    The Ohio State University. Columbus, OH: 2024.
#' @source For more details, please visit the offical website of the dataset \url{https://www.nlsinfo.org/}
#' @examples
#'  res.mari <- panelTVP(Used.Mari.Since.DLI ~ Sex + Baseline.Age_c + Bullied.Before.12 + Ethnicity + Residence,
#'                       data = Marijuana,
#'                       t = Marijuana$t,
#'                       id = Marijuana$id,
#'                       mcmc.opt = list(chain.length = 100, burnin = 0, thin = 1, asis = TRUE),
#'                       model = "Logit")
#'  # Please increase the length of the Markov Chain. 100 draws in total is only for demonstration!
#'
"Marijuana"

#' Income dataset for Gaussian panel models
#'
#' This dataset is a pre-processed subset of the National Longitudinal Survey of Youth 1997,
#' which is sponsored by the U.S. Bureau of Labor Statistics. This panel study provides a
#' rich data source on a wide range of socio-economic topics. This subsetted dataset serves
#' mainly as an application for Gaussian time-varying panel models. The goal of this
#' dataset is to analyze the gross family income (on the log-scale) in terms of covariates.
#' Note that this is a subsetted version of the official dataset that only contains data
#' of 2015 subjects form surveys conducted in 2015, 2017, 2019 and 2021.
#' The actual study began in 1997 and
#' comprises a total of 20 waves.
#'
#' \describe{
#'   \item{Sex}{biological sex of the respondent (factor variable)}
#'   \item{Ethnicity}{ethnicity of the respondent (factor variable)}
#'   \item{Year}{year of survey round (numeric variable)}
#'   \item{Edu}{the person's highest completed level of formal education at interview date}
#'   \item{Hours.Worked}{indicator whether the person has worked 30 hours or more
#'     per week (factor variable)}
#'   \item{Age}{age of respondent (numeric variable)}
#'   \item{t}{time index (numeric variable)}
#'   \item{Gross.Family.Income_log}{gross family income of the last year on
#'   the log-scale (numeric variable)}
#'   \item{id}{respondent identifier (numeric variable); note that some id's
#'    are missing as some original respondents did not met the inclusion
#'    requirements}
#'   \item{Baseline.Age_c}{centered baseline age, i.e., a person with a value of 0
#'    was 30 years old in 2015 (numeric variable)}
#' }
#' @references Bureau of Labor Statistics, U.S. Department of Labor.
#'  National Longitudinal Survey of Youth 1997 cohort, 1997-2021 (rounds 1-20).
#'   Produced and distributed by the Center for Human Resource Research (CHRR),
#'    The Ohio State University. Columbus, OH: 2024.
#' @source For more details, please visit the offical website of the dataset \url{https://www.nlsinfo.org/}
#' @examples
#'  res.income <- panelTVP(Gross.Family.Income_log ~ Sex + Baseline.Age_c + Ethnicity + Edu + Hours.Worked,
#'                       data = Income,
#'                       t = Income$t,
#'                       id = Income$id,
#'                       mcmc.opt = list(chain.length = 100, burnin = 0, thin = 1, asis = TRUE),
#'                       model = "Gaussian")
#'  # Please increase the length of the Markov Chain. 100 draws in total is only for demonstration!
#'
"Income"
