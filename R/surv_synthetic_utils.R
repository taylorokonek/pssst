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