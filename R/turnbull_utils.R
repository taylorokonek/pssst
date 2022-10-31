#' Function to input to optim, calls rcpp_loglik_multi
#' 
#' @param data dataframe
#' @param weights vector of survey weights
#' @param par parameter vector
#' @param shape_par_ids which pars are shape parameters
#' @param dist integer denoting which distribution to use
#' @return negative log likelihood
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
