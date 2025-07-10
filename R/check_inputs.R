# check input arguments of sim_panelTVP ----------------------------------------

check_sim <- function(n,
                      Tmax,
                      model,
                      beta = NULL,
                      theta = NULL,
                      lambda = NULL,
                      psi = NULL,
                      r = NULL,
                      sigma2 = NULL,
                      beta.nb = NULL,
                      theta.nb = NULL,
                      lambda.nb = NULL,
                      psi.nb = NULL,
                      beta.logit = NULL,
                      theta.logit = NULL,
                      lambda.logit = NULL,
                      psi.logit = NULL){

  if(!is.numeric(n) || length(n) != 1 || n %% 1 != 0 || !is.finite(n)){
    stop("Argument 'n' must be a single value that represents the number of subjects")
  }
  if(!is.numeric(Tmax) || length(Tmax) != 1 || Tmax %% 1 != 0 || !is.finite(Tmax) || Tmax <= 2){
    stop("Argument 'Tmax' must be a single value > 2 that represents the number of repeated measurements")
  }
  if(!(model %in% c("Gaussian", "Probit", "Logit", "NegBin", "ZINB"))){
    stop("Argument 'model' must be either 'Gaussian', 'Probit', 'Logit', 'NegBin' or 'ZINB'")
  }
  if(model %in% c("NegBin", "ZINB") && is.null(r)){
    stop("Argument 'r' must be specified when model is either 'NegBin' or 'ZINB'")
  }
  if(model == "Gaussian" && is.null(sigma2)){
    stop("Argument 'sigma2' must be specified when model is 'Gaussian'")
  }
  if(model %in% c("Gaussian", "Probit", "Logit", "NegBin") &&
     any(sapply(list(beta, theta, lambda, psi), is.null))){
    stop("When model is either 'Gaussian', 'Probit', 'Logit', or 'NegBin' arguments 'beta', 'theta', 'lambda' and 'psi' have to be specified")
  }
  if(model == "ZINB" && any(sapply(list(beta.nb, theta.nb, lambda.nb, psi.nb,
                                        beta.logit, theta.logit, lambda.logit, psi.logit), is.null))){
    stop("When model is 'ZINB' arguments 'beta.nb', 'theta.nb', 'lambda.nb', 'psi.nb', 'beta.logit', 'theta.logit', 'lambda.logit', 'psi.logit' have to be specified")
  }
  if(all(sapply(list(beta, theta, lambda, psi), function(x) !is.null(x)))){
    if(any(sapply(list(beta, theta, lambda, psi), function(x) length(x) == 0))){
      stop("Arguments 'beta', 'theta', 'lambda', 'psi' must be of positive length")
    }
    if(length(beta) != length(theta)) stop("Argument 'beta' must be of same length as argument 'theta'")
    if(!is.numeric(beta) || any(!is.finite(beta))) stop("Argument 'beta' must be vector of finite values")
    if(!is.numeric(theta) || any(!is.finite(theta))) stop("Argument 'theta' must be vector of finite values")
    if(!is.numeric(lambda) || length(lambda) > 1 || !is.finite(lambda)) stop("Argument 'lambda' must be finite numeric")
    if(!is.numeric(psi) || length(psi) > 1 || !is.finite(psi)) stop("Argument 'psi' must be finite numeric")
  }
  if(all(sapply(list(beta.nb, theta.nb, lambda.nb, psi.nb,
                     beta.logit, theta.logit, lambda.logit, psi.logit), function(x) !is.null(x)))){
    if(any(sapply(list(beta.nb, theta.nb, lambda.nb, psi.nb,
                       beta.logit, theta.logit, lambda.logit, psi.logit), function(x) length(x) == 0))){
      stop("Arguments 'beta.nb', 'theta.nb', 'lambda.nb', 'psi.nb', 'beta.logit', 'theta.logit', 'lambda.logit', 'psi.logit' must be of positive length")
    }
    if(length(beta.nb) != length(theta.nb)) stop("Argument 'beta.nb' must be of same length as argument 'theta.nb'")
    if(!is.numeric(beta.nb) || any(!is.finite(beta.nb))) stop("Argument 'beta.nb' must be vector of finite values")
    if(!is.numeric(theta.nb) || any(!is.finite(theta.nb))) stop("Argument 'theta.nb' must be vector of finite values")
    if(!is.numeric(lambda.nb) || length(lambda.nb) > 1 || !is.finite(lambda.nb)) stop("Argument 'lambda.nb' must be finite numeric")
    if(!is.numeric(psi.nb) || length(psi.nb) > 1 || !is.finite(psi.nb)) stop("Argument 'psi.nb' must be finite numeric")
    if(length(beta.logit) != length(theta.logit)) stop("Argument 'beta.logit' must be of same length as argument 'theta.logit'")
    if(!is.numeric(beta.logit) || any(!is.finite(beta.logit))) stop("Argument 'beta.logit' must be vector of finite values")
    if(!is.numeric(theta.logit) || any(!is.finite(theta.logit))) stop("Argument 'theta.logit' must be vector of finite values")
    if(!is.numeric(lambda.logit) || length(lambda.logit) > 1 || !is.finite(lambda.logit)) stop("Argument 'lambda.logit' must be finite numeric")
    if(!is.numeric(psi.logit) || length(psi.logit) > 1 || !is.finite(psi.logit)) stop("Argument 'psi.logit' must be finite numeric")
  }
  if(!is.null(r) && (!is.numeric(r) || length(r) > 1 || r < 0 || !is.finite(r))){
    stop("Argument 'r' must be a finite, positive numeric")
  }
  if(!is.null(sigma2) && (!is.numeric(sigma2) || length(sigma2) > 1 || sigma2 < 0 || !is.finite(sigma2))){
    stop("Argument 'sigma2' must be a finite, positive numeric")
  }

}

# check input arguments of panelTVP --------------------------------------------

# check input arguments of predict S3-functions --------------------------------

check.predict <- function(model, X.new, timepoint, coverage, pop.pred){

  if(!(is.matrix(X.new) || is.data.frame(X.new))){
    stop("Argument 'X.new' must be either a matrix or a data frame.")
  }

  if(ncol(X.new) != ncol(model$data$X)){
    stop("Number of columns in X.new must match number of columns in original design matrix. In case you fitted an intercept, X.new must contain a column of ones.")
  }

  if(is.null(colnames(X.new))){
    warning("Columns of argument 'X.new' are not labelled. This is dangerous as you must make sure that the columns match the columns of the original design matrix.")
  }

  else{

    if(any(!(colnames(X.new) %in% colnames(model$data$X)))){
      stop("Variable names of X.new must match the corresponding names from the original design matrix.")
    }

    if(sum(colnames(X.new) != colnames(model$data$X)) > 0){
      warning("Column names in 'X.new' were reordered to match the corresponding ones from the original design matrix.")
      X.new <- X.new[, colnames(model$data$X), drop = FALSE] # reordering variables w.r.t. original design matrix
    }

  }

  if(!is.numeric(timepoint) || length(timepoint) != 1 || timepoint %% 1 != 0 || !is.finite(timepoint)){
    stop("Argument 'timepoint' must be a single numeric value. If you want predictions for multiple time points, please call the predict function multiple times.")
  }

  if(!is.numeric(coverage) || length(coverage) != 1 || !is.finite(coverage) || coverage < 0 || coverage > 1){
    stop("Argument 'coverage' must be a single numeric value between 0 and 1.")
  }

  if(!is.logical(pop.pred) || length(pop.pred) != 1 || !is.finite(pop.pred)){
    stop("Argument 'pop.pred' must be a single logical value.")
  }

}

check.predict_ZINB <- function(model, X_nb.new, X_logit.new, timepoint, coverage, pop.pred){

  if(!(is.matrix(X_nb.new) || is.data.frame(X_nb.new))){
    stop("Argument 'X_nb.new' must be either a matrix or a data frame.")
  }

  if(ncol(X_nb.new) != ncol(model$data$X_nb)){
    stop("Number of columns in X_nb.new must match number of columns in original design matrix. In case you fitted an intercept, X_nb.new must contain a column of ones.")
  }

  if(is.null(colnames(X_nb.new))){
    warning("Columns of argument 'X_nb.new' are not labelled. This is dangerous as you must make sure that the columns match the columns of the original design matrix.")
  }

  else{

    if(any(!(colnames(X_nb.new) %in% colnames(model$data$X_nb)))){
      stop("Variable names of X_nb.new must match the corresponding names from the original design matrix.")
    }

    if(sum(colnames(X_nb.new) != colnames(model$data$X_nb)) > 0){
      warning("Column names in 'X_nb.new' were reordered to match the corresponding ones from the original design matrix.")
      X_nb.new <- X_nb.new[, colnames(model$data$X_nb), drop = FALSE] # reordering variables w.r.t. original design matrix
    }

  }

  if(!(is.matrix(X_logit.new) || is.data.frame(X_logit.new))){
    stop("Argument 'X_logit.new' must be either a matrix or a data frame.")
  }

  if(ncol(X_logit.new) != ncol(model$data$X_logit)){
    stop("Number of columns in X_logit.new must match number of columns in original design matrix. In case you fitted an intercept, X_logit.new must contain a column of ones.")
  }

  if(is.null(colnames(X_logit.new))){
    warning("Columns of argument 'X_logit.new' are not labelled. This is dangerous as you must make sure that the columns match the columns of the original design matrix.")
  }

  else{

    if(any(!(colnames(X_logit.new) %in% colnames(model$data$X_logit)))){
      stop("Variable names of X_logit.new must match the corresponding names from the original design matrix.")
    }

    if(sum(colnames(X_logit.new) != colnames(model$data$X_logit)) > 0){
      warning("Column names in 'X_logit.new' were reordered to match the corresponding ones from the original design matrix.")
      X_logit.new <- X_logit.new[, colnames(model$data$X_logit), drop = FALSE] # reordering variables w.r.t. original design matrix
    }

  }

  if(!is.numeric(timepoint) || length(timepoint) != 1 || timepoint %% 1 != 0 || !is.finite(timepoint)){
    stop("Argument 'timepoint' must be a single numeric value. If you want predictions for multiple time points, please call the predict function multiple times.")
  }

  if(!is.numeric(coverage) || length(coverage) != 1 || !is.finite(coverage) || coverage < 0 || coverage > 1){
    stop("Argument 'coverage' must be a single numeric value between 0 and 1.")
  }

  if(!is.logical(pop.pred) || length(pop.pred) != 1 || !is.finite(pop.pred)){
    stop("Argument 'pop.pred' must be a single logical value.")
  }

}
