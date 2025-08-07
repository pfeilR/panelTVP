#' @title Get plots of time-varying parameters based on a \code{panelTVP.Gaussian} object
#'
#' @description
#'  This \code{plot} function produces effect plots for the time-varying parameters,
#'  i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T} separately for each covariate. Plots are produced
#'   using the \code{ggplot2} library.
#'
#' @param x an object of class \code{panelTVP.Gaussian}
#' @param nplots indicates how many plots should be printed on one page
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method plot panelTVP.Gaussian
#' @references
#'  Wickham, H. (2016). \code{ggplot2}: Elegant Graphics for Data Analysis.
#'  Springer Verlag, New York.
#' @examples
#' # Plot results from an object of class panelTVP.Gaussian
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.gaussian <- sim_panelTVP(n = 100,
#'                              Tmax = 4,
#'                              beta = c(4,1,0,0),
#'                              theta = c(1,0.5,0,0),
#'                              lambda = 1,
#'                              psi = 0.2,
#'                              model = "Gaussian",
#'                              sigma2 = 0.7)
#' res.gaussian <- panelTVP(y ~ W1 + W2 + W3,
#'                          data = sim.gaussian$observed,
#'                          id = sim.gaussian$observed$id,
#'                          t = sim.gaussian$observed$t,
#'                          mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                          model = "Gaussian")
#' plot(res.gaussian, nplots = 1)
plot.panelTVP.Gaussian <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters for a \code{panelTVP.Probit} object
#'
#' @description
#'  This \code{plot} function produces effect plots for the time-varying parameters,
#'  i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T} separately for each covariate. Plots are produced
#'   using the \code{ggplot2} library.
#'
#' @param x an object of class \code{panelTVP.Probit}
#' @param nplots indicates how many plots should be printed on one page
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method plot panelTVP.Probit
#' @references
#'  Wickham, H. (2016). \code{ggplot2}: Elegant Graphics for Data Analysis.
#'  Springer Verlag, New York.
#' @examples
#' # Plot results from an object of class panelTVP.Probit
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.probit <- sim_panelTVP(n = 100,
#'                            Tmax = 4,
#'                            beta = c(1,0.5,0,0),
#'                            theta = c(0.8,0.5,0,0),
#'                            lambda = 1,
#'                            psi = 0.2,
#'                            model = "Probit")
#' res.probit <- panelTVP(y ~ W1 + W2 + W3,
#'                        data = sim.probit$observed,
#'                        id = sim.probit$observed$id,
#'                        t = sim.probit$observed$t,
#'                        mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                        model = "Probit")
#' plot(res.probit, nplots = 1)
plot.panelTVP.Probit <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters based on a \code{panelTVP.Logit} object
#'
#' @description
#'  This \code{plot} function produces effect plots for the time-varying parameters,
#'  i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T} separately for each covariate. Plots are produced
#'   using the \code{ggplot2} library.
#'
#' @param x an object of class \code{panelTVP.Logit}
#' @param nplots indicates how many plots should be printed on one page
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method plot panelTVP.Logit
#' @references
#'  Wickham, H. (2016). \code{ggplot2}: Elegant Graphics for Data Analysis.
#'  Springer Verlag, New York.
#' @examples
#' # Plot results from an object of class panelTVP.Logit
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.logit <- sim_panelTVP(n = 100,
#'                           Tmax = 4,
#'                           beta = c(1,0.5,0,0),
#'                           theta = c(0.8,0.5,0,0),
#'                           lambda = 1,
#'                           psi = 0.2,
#'                           model = "Logit")
#' res.logit <- panelTVP(y ~ W1 + W2 + W3,
#'                       data = sim.logit$observed,
#'                       id = sim.logit$observed$id,
#'                       t = sim.logit$observed$t,
#'                       mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                       model = "Logit")
#' plot(res.logit, nplots = 1)
plot.panelTVP.Logit <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters based on a \code{panelTVP.NegBin} object
#'
#' @description
#'  This \code{plot} function produces effect plots for the time-varying parameters,
#'  i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T} separately for each covariate. Plots are produced
#'   using the \code{ggplot2} library.
#'
#' @param x an object of class \code{panelTVP.NegBin}
#' @param nplots indicates how many plots should be printed on one page
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method plot panelTVP.NegBin
#' @references
#'  Wickham, H. (2016). \code{ggplot2}: Elegant Graphics for Data Analysis.
#'  Springer Verlag, New York.
#' @examples
#' # Plot results from an object of class panelTVP.NegBin
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.negbin <- sim_panelTVP(n = 100,
#'                            Tmax = 4,
#'                            beta = c(1,0.5,0,0),
#'                            theta = c(0.8,0.5,0,0),
#'                            lambda = 1,
#'                            psi = 0.2,
#'                            r = 2,
#'                            model = "NegBin")
#' res.negbin <- panelTVP(y ~ W1 + W2 + W3,
#'                        data = sim.negbin$observed,
#'                        id = sim.negbin$observed$id,
#'                        t = sim.negbin$observed$t,
#'                        mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                        model = "NegBin")
#' plot(res.negbin, nplots = 1)
plot.panelTVP.NegBin <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters based on a \code{panelTVP.ZINB} object
#'
#' @description
#'  This \code{plot} function produces effect plots for the time-varying parameters,
#'  i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T} separately for each covariate. Plots are produced
#'   using the \code{ggplot2} library. Note that this
#'   function either plots the effects of the Negative Binomial (or count) model, or
#'   the Logit (or zero-inflation) model.
#'
#' @param x an object of class \code{panelTVP.ZINB}
#' @param component either 'count' or 'inflation' to create plots for either the count or the
#'  zero-inflation component (no default)
#' @param nplots indicates how many plots should be printed on one page
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler, Helga Wagner
#' @exportS3Method plot panelTVP.ZINB
#' @references
#'  Wickham, H. (2016). \code{ggplot2}: Elegant Graphics for Data Analysis.
#'  Springer Verlag, New York.
#' @examples
#' # Plot results from an object of class panelTVP.ZINB
#' # NB: To reduce computational effort, we have drastically reduced the length
#' # of the Markov Chain. You should use a much longer chain in your applications.
#' sim.zinb <- sim_panelTVP(n = 100,
#'                          Tmax = 4,
#'                          beta_zinb.count = c(0.5,-0.7,0,0),
#'                          theta_zinb.count = c(0.05,0.5,0,0),
#'                          lambda_zinb.count = 0.5,
#'                          psi_zinb.count = 0.02,
#'                          beta_zinb.inflation = c(-1,0.6,0,0),
#'                          theta_zinb.inflation = c(0,1,0,0),
#'                          lambda_zinb.inflation = 0.7,
#'                          psi_zinb.inflation = 0,
#'                          r = 2,
#'                          model = "ZINB")
#' res.zinb <- panelTVP(y ~ W1.nb + W2.nb + W3.nb | W1.logit + W2.logit + W3.logit,
#'                      data = sim.zinb$observed,
#'                      id = sim.zinb$observed$id,
#'                      t = sim.zinb$observed$t,
#'                      mcmc.opt = list(chain.length = 200, burnin = 100, thin = 1, asis = TRUE),
#'                      model = "ZINB")
#' plot(res.zinb, nplots = 1, component = "count") # count component
#' plot(res.zinb, nplots = 1, component = "inflation") # zero-inflation component
plot.panelTVP.ZINB <- function(x, component = NULL, nplots = 4, ...){
  if(is.null(component) || length(component) != 1 || !(component %in% c("count", "inflation"))){
    stop("Argument 'component' must be either 'count' or 'inflation' depending on which parameters you are interested in.")
  }
  if(component == "count"){
    return(plot_effects(summary_table = x$posterior_nb, Tmax = x$data$Tmax, X = x$data$X_nb, nplots = nplots))
  } else{
    return(plot_effects(summary_table = x$posterior_logit, Tmax = x$data$Tmax, X = x$data$X_logit, nplots = nplots))
  }
}


plot_effects <- function(summary_table, Tmax, X, nplots = 4){

  beta_lambda_rows <- c(grep("^beta_t\\d+", rownames(summary_table), value = TRUE),
                        grep("^lambda_t[0-9]+", rownames(summary_table), value = TRUE))
  beta_lambda_data <- summary_table[beta_lambda_rows, , drop = FALSE]
  if(sum(startsWith(rownames(summary_table), "lambda_t")) == 1){ # cps
    beta_lambda_rows <- c(beta_lambda_rows, paste0("lambda_t",1:Tmax))
    mat <- matrix(summary_table[startsWith(rownames(summary_table), "lambda_t"),], nrow = 1)
    mat <- mat[rep(1,Tmax),]
    rownames(mat) <- paste0("lambda_t", 1:Tmax)
    beta_lambda_data <- rbind(beta_lambda_data, mat)
  }
  index_df <- matrix(nrow = (ncol(X)+1)*Tmax, ncol = 2)
  index_df[,1] <- c(rep(paste0("x",1:ncol(X)), each = Tmax), rep("Factor Loading", Tmax))
  index_df[,2] <- rep(1:Tmax, ncol(X)+1)
  index_df <- as.data.frame(index_df)
  colnames(index_df) <- c("var", "time")
  index_df$time <- as.numeric(index_df$time)
  beta_lambda_df <- data.frame(
    name = beta_lambda_rows,
    lower = beta_lambda_data[, "LO"],
    mean  = beta_lambda_data[, "mean"],
    upper = beta_lambda_data[, "UP"]
  )
  true_names <- c(colnames(X), "Random Effect")
  index_df$varname <- rep(true_names, each = Tmax)
  plot_df <- cbind(index_df, beta_lambda_df)
  plot_list <- dplyr::group_split(dplyr::group_by(plot_df, var))
  plot_list <- plot_list[order(
    sapply(plot_list, function(x){
      var_value <- x$var[1]
      if(var_value == "Factor Loading"){
        return(0)
      }else{
        return(as.numeric(sub("x", "", var_value)))
      }
    })
  )]
  plot_objs <- list()
  for(i in 1:length(plot_list)) plot_objs[[i]] <- make_plot(plot_list[[i]])
  plot_objs[[1]] <- plot_objs[[1]] + ggplot2::ylab(expression(hat(lambda)))
  n_total <- length(plot_objs)
  i <- 1
  while(i <= n_total) {
    end <- min(i + nplots - 1, n_total)
    grid_plots <- plot_objs[i:end]
    gridExtra::grid.arrange(grobs = grid_plots, ncol = 1)
    i <- end + 1
    if(i <= n_total) readline("The plots are out there, hit [Return] to see ...")
  }

}

utils::globalVariables(c("lower", "upper", "mean", "time", "varname"))

make_plot <- function(df){
  range_y <- range(df$lower, df$upper, na.rm = TRUE)
  y_pad <- 0.1 * diff(range_y)
  y_min <- min(df$lower, na.rm = TRUE) - y_pad
  y_max <- max(df$upper, na.rm = TRUE) + y_pad
  ggplot2::ggplot(df, ggplot2::aes(x = time, y = mean)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "#A6CEE3", alpha = 0.5) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_point(color = "black", size = 2) +
    ggplot2::scale_x_continuous(breaks = df$time) +
    ggplot2::scale_y_continuous(limits = c(y_min, y_max)) +
    ggplot2::labs(
      title = unique(df$varname[1]),
      x = "Timepoint", y = expression(hat(beta))
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major = ggplot2::element_line(color = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      text = ggplot2::element_text(color = "black")
    )
}
