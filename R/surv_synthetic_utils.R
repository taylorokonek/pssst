#' Function to input to optim, calls rcpp_loglik_multi
#' 
#' @param par parameter vector
#' @param data dataframe
#' @param weights vector of survey weights
#' @param shape_par_ids which pars are shape parameters
#' @param dist integer denoting which distribution to use
#' @param num_periods number of periods in dataframe
#' @param etsp_c if dist is etsp, fixed value to use for the c parameter. Default is zero.
#' @return negative log likelihood
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
optim_fn <- function(par, data, weights, shape_par_ids, dist, breakpoints,
                     num_periods, etsp_c = 0) {
  # weibull
  if (dist == 1) {
    pars_per_period <- length(par[-shape_par_ids]) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)

        a <- rcpp_loglik_multi(x_df = data, 
                           num_periods = num_periods,
                           log_shapes = par[shape_par_ids], 
                           log_scales = par[-shape_par_ids], 
                           dist = dist,
                           breakpoints = breakpoints,
                           par_period_id = par_period_id)
    
    ret <- -sum(a * weights)
    
    # exponential, piecewise exponential
  } else if (dist == 0) {
    pars_per_period <- length(par) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    # can fit exponential likelihood using weibull code with shapes == 1
    ret <- -sum(rcpp_loglik_multi(x_df = data, 
                                  num_periods = num_periods,
                                  log_shapes = 0, 
                                  log_scales = par, 
                                  dist = 1,
                                  breakpoints = breakpoints,
                                  par_period_id = par_period_id) * weights)
    
  } else if (dist == 2) {
    pars_per_period <- length(par) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    
    a <- rcpp_loglik_multi(x_df = data, 
                           num_periods = num_periods,
                           log_shapes = 0, 
                           log_scales = par, 
                           dist = dist,
                           breakpoints = breakpoints,
                           par_period_id = par_period_id)
    ret <- -sum(a * weights)
  } else if (dist == 3) {
    pars_per_period <- length(par[-shape_par_ids]) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    a <- rcpp_loglik_multi(x_df = data, 
                                  num_periods = num_periods,
                                  log_shapes = par[shape_par_ids], 
                                  log_scales = par[-shape_par_ids], 
                                  dist = dist,
                                  breakpoints = breakpoints,
                                  par_period_id = par_period_id) 
    ret <- -sum(a * weights)
    
  } else if (dist == 4 | dist == 5 | dist == 6 | dist == 7 | dist == 8 ) {
    pars_per_period <- length(par[-shape_par_ids]) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)

    a <- rcpp_loglik_multi(x_df = data, 
                           num_periods = num_periods,
                           log_shapes = par[shape_par_ids], 
                           log_scales = par[-shape_par_ids], 
                           dist = dist,
                           breakpoints = breakpoints,
                           par_period_id = par_period_id,
                           etsp_c = etsp_c) 
    ret <- -sum(a * weights)

  } 
  return(ret)
}

#' Function to input to grad, calls rcpp_loglik_multi
#' 
#' @param par parameter vector
#' @param data dataframe
#' @param weights vector of survey weights
#' @param shape_par_ids which pars are shape parameters
#' @param dist integer denoting which distribution to use
#' @param num_periods numer of periods in dataframe
#' @param etsp_c if dist is etsp, fixed value to use for the c parameter. Default is zero.
#' @return negative log likelihood
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
optim_fn_grad <- function(par, data, weights, shape_par_ids, dist,
                          num_periods, breakpoints, etsp_c = 0) {
  if (dist == 1) {
    pars_per_period <- length(par[-shape_par_ids]) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    a <- sum(rcpp_loglik_multi(x_df = data, 
                          num_periods = num_periods,
                          log_shapes = par[shape_par_ids], 
                          log_scales = par[-shape_par_ids], 
                          dist = dist,
                          breakpoints = breakpoints,
                          par_period_id = par_period_id) * weights)
  } else if (dist == 0) {
    pars_per_period <- length(par) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    # can fit exponential likelihood using weibull code with shapes == 1
    a <- sum(rcpp_loglik_multi(x_df = data, 
                          num_periods = num_periods,
                          log_shapes = 0, 
                          log_scales = par, 
                          dist = 1,
                          breakpoints = breakpoints,
                          par_period_id = par_period_id) * weights)
  } else if (dist == 2) {
    pars_per_period <- length(par) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    a <- sum(rcpp_loglik_multi(x_df = data, 
                          num_periods = num_periods,
                          log_shapes = 0, 
                          log_scales = par, 
                          dist = dist,
                          breakpoints = breakpoints,
                          par_period_id = par_period_id) * weights)
  } else if (dist == 3) {
    pars_per_period <- length(par[-shape_par_ids]) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    a <- sum(rcpp_loglik_multi(x_df = data, 
                           num_periods = num_periods,
                           log_shapes = par[shape_par_ids], 
                           log_scales = par[-shape_par_ids], 
                           dist = dist,
                           breakpoints = breakpoints,
                           par_period_id = par_period_id) * weights)
  } else if (dist == 4 | dist == 5 | dist == 6 | dist == 7 | dist == 8) {
    pars_per_period <- length(par[-shape_par_ids]) / num_periods
    par_period_id <- rep(1:num_periods, each = pars_per_period)
    
    a <- sum(rcpp_loglik_multi(x_df = data, 
                               num_periods = num_periods,
                               log_shapes = par[shape_par_ids], 
                               log_scales = par[-shape_par_ids], 
                               dist = dist,
                               breakpoints = breakpoints,
                               par_period_id = par_period_id,
                               etsp_c = etsp_c) * weights)
    
  }
  return(a)
}

#' CDF for piecewise exponential distribution
#'
#' @param x value
#' @param log_scales log scales
#' @param breakpoints breakpoints
#' @param log return log cdf or no
#' @return cdf
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
p_piecewise_exponential <- function(x, log_scales, breakpoints, log = FALSE) {
  # F(x) = 1 - exp(-H(x))
  H_x <- rcpp_hazard_integral(0, x, 1, log_scales, 2, breakpoints)
  if (log) {
    log(1 - exp(-H_x))
  } else {
    1 - exp(-H_x)
  }
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

#' derivative of CDF for exponential distribution wrt log_scale
#' 
#' @param x value
#' @param log_scale log_scale
#' @return derivative wrt log_scale of exponential cdf 
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
d_p_exponential_log_scale <- function(x, log_scale) {
  -(exp(-(x/exp(log_scale))) * (x * exp(log_scale)/exp(log_scale)^2))
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

#' cdf delta method for exponential distribution 
#' can be used for u5mr (x = 60), imr (x = 12), nmr (x = 1), etc.
#' 
#' @param x value
#' @param par c(log_scale)
#' @param vmat vcov matrix for par
#' @return delta method variance for cdf evaluated at x
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
p_exponential_deltamethod <- function(x, par, vmat) {
  grad <- diag(length(par)) * d_p_exponential_log_scale(x, par)
  delta_vmat <- t(grad) %*% vmat %*% grad
  return(delta_vmat)
}

#' cdf delta method for piecewise exponential distribution 
#' can be used for u5mr (x = 60), imr (x = 12), nmr (x = 1), etc.
#' 
#' @param x value
#' @param par c(log_scale)
#' @param vmat vcov matrix for par
#' @return delta method variance for cdf evaluated at x
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
p_piecewise_exponential_deltamethod <- function(x, par, vmat) {
  grad <- diag(length(par)) * d_p_exponential_log_scale(x, par)
  delta_vmat <- t(grad) %*% vmat %*% grad
  return(delta_vmat)
}



#' h(x) for ETSP
#' 
#' @param x value
#' @param a parameter
#' @param b parameter
#' @param c parameter
#' @param p parameter
#' @return hazard for ETSP
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
h_etsp <- function (x, a, b, c, p) {
  a * (x + c)^(-p) * exp(-b * x)
}



#' H(x) for ETSP
#' 
#' @param x value
#' @param lower lower val for integral of hazard
#' @param upper upper val for integral of hazard
#' @param a upper val for integral of hazard
#' @param b upper val for integral of hazard
#' @param c upper val for integral of hazard
#' @param p upper val for integral of hazard
#' @return cumulative hazard for ETSP family of hazards between lower and upper 
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
H_etsp <- function(lower, upper, a, b, c, p) {
  ret <- -1*a*exp(b*c)*b^(p-1)*(expint::gammainc(1-p, b*(c+upper)) - expint::gammainc(1-p, b*(c+lower)))
  return(ret)
}








