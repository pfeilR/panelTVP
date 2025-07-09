#' @title Get plots of time-varying parameters for an object of class \code{panelTVP.Gaussian}
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
#' @author Roman Pfeiler
#' @exportS3Method plot panelTVP.Gaussian
plot.panelTVP.Gaussian <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters for an object of class \code{panelTVP.Probit}
#'
#' @description
#'  This \code{plot} function produces effect plots for the time-varying parameters,
#'  i.e., \eqn{\boldsymbol{\beta}_1,\dots,\boldsymbol{\beta}_T} and
#'   \eqn{\lambda_1,\dots,\lambda_T} separately for each covariate Plots are produced
#'   using the \code{ggplot2} library.
#'
#' @param x an object of class \code{panelTVP.Probit}
#' @param nplots indicates how many plots should be printed on one page
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @exportS3Method plot panelTVP.Probit
plot.panelTVP.Probit <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters for an object of class \code{panelTVP.Logit}
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
#' @author Roman Pfeiler
#' @exportS3Method plot panelTVP.Logit
plot.panelTVP.Logit <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters for an object of class \code{panelTVP.NegBin}
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
#' @author Roman Pfeiler
#' @exportS3Method plot panelTVP.NegBin
plot.panelTVP.NegBin <- function(x, nplots = 4, ...){
  plot_effects(summary_table = x$posterior, Tmax = x$data$Tmax, X = x$data$X, nplots = nplots)
}

#' @title Get plots of time-varying parameters for an object of class \code{panelTVP.ZINB}
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
#' @param component either 'NegBin' or 'Logit' to create plots for either the count or the
#'  zero-inflation component, respectively - must be specified!
#' @param nplots indicates how many plots should be printed on one page
#' @param ... optional arguments passed to the function (those are ignored)
#'
#' @author Roman Pfeiler
#' @exportS3Method plot panelTVP.ZINB
plot.panelTVP.ZINB <- function(x, component = NULL, nplots = 4, ...){
  if(is.null(component) || length(component) != 1 || !(component %in% c("NegBin", "Logit"))){
    stop("Argument 'component' must be either NegBin or Logit depending on which parameters you are interested in.")
  }
  if(component == "NegBin"){
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
  parse_index <- function(name){
    if(grepl("^beta_t\\d{2,}$", name)){
      matches <- regmatches(name, regexec("^beta_t(\\d)(\\d+)$", name))[[1]]
      return(data.frame(
        var = paste0("x", matches[2]),
        time = as.integer(matches[3]),
        stringsAsFactors = FALSE
      ))
    }
    if(grepl("^lambda_t\\d+$", name)){
      matches <- regmatches(name, regexec("^lambda_t(\\d+)$", name))[[1]]
      return(data.frame(
        var = "Factor Loading",
        time = as.integer(matches[2]),
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  }
  index_df <- do.call(rbind, lapply(beta_lambda_rows, parse_index))
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
