#' Function to input to optim, calls rcpp_loglik_multi
#' 
#' @param data dataframe
#' @param weights vector of survey weights
#' @param par parameter vector
#' @param shape_par_ids which pars are shape parameters
#' @return negative log likelihood
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
optim_fn <- function(data, weights, par, shape_par_ids) {
  -sum(rcpp_loglik_multi(data, par[shape_par_ids], par[-shape_par_ids]) * weights)
}

#' CDF for weibull distribution
#' 
#' @param x value
#' @param log_shape log_shape
#' @param log_scale log_scale
#' @param log return log cdf or no
#' @return cdf
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
p_weibull <- function(x, log_shape, log_scale, log = FALSE) {
  if (log) {
    log(1 - exp(-(x/exp(log_scale))^(exp(log_shape))))
  } else {
    1 - exp(-(x/exp(log_scale))^(exp(log_shape)))
  }
}

#' derivative of CDF for weibull distribution wrt log_scale
#' 
#' @param x value
#' @param log_shape log_shape
#' @param log_scale log_scale
#' @return derivative wrt log_scale of weibull cdf 
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
d_p_weibull_log_scale <- function(x, log_shape, log_scale) {
  -(exp(-(x/exp(log_scale))^(exp(log_shape))) * ((x/exp(log_scale))^((exp(log_shape)) - 
                                                                       1) * ((exp(log_shape)) * (x * exp(log_scale)/exp(log_scale)^2))))
}

#' derivative of CDF for weibull distribution wrt log_shape
#' 
#' @param x value
#' @param log_shape log_shape
#' @param log_scale log_scale
#' @return derivative wrt log_shape of weibull cdf 
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
d_p_weibull_log_shape <- function(x, log_shape, log_scale) {
  exp(-(x/exp(log_scale))^(exp(log_shape))) * ((x/exp(log_scale))^(exp(log_shape)) * 
                                                 (log((x/exp(log_scale))) * exp(log_shape)))
}

#' cdf delta method for weibull distribution 
#' can be used for u5mr (x = 60), imr (x = 12), nmr (x = 1), etc.
#' 
#' @param x value
#' @param par c(log_shape, log_scale)
#' @param vmat vcov matrix for par
#' @param shape_par_ids which parameters are shapes
#' @return delta method variance for cdf evaluated at x
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
p_weibull_deltamethod <- function(x, par, vmat, shape_par_ids) {
  shapes <- par[shape_par_ids]
  scales <- par[-shape_par_ids]
  
  if (length(shape_par_ids) == 1) {
    g_log_shape <- d_p_weibull_log_shape(x, shapes, scales)
    g_log_scale <- diag(length(scales)) * d_p_weibull_log_scale(x, shapes, scales)
  } else {
    g_log_shape <- diag(length(shapes)) * d_p_weibull_log_shape(x, shapes, scales)
    g_log_scale <- diag(length(scales)) * d_p_weibull_log_scale(x, shapes, scales)
  }
  
  grad <- rbind(g_log_shape, g_log_scale)
  delta_vmat <- t(grad) %*% vmat %*% grad
  
  return(delta_vmat)
}








