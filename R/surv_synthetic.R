#' Fit a parametric survival model to DHS data across multiple time periods.
#' 
#' Fits a parametric survival model to DHS data (formatted from \code{format_dhs}) to get survey-weighted
#' parameter estimates and finite population variances for parametric survival curves in each time period.
#' 
#' @param df a dataframe containing the output from \code{format_dhs}, or optionally, dataframe
#' containing the following columns
#' @param individual column corresponding to individual ID in \code{df}
#' @param household column corresponding to household ID in \code{df}
#' @param cluster column corresponding to cluster ID in \code{df}
#' @param strata column corresponding to strata ID in \code{df}
#' @param weights column corresponding to weights in \code{df}
#' @param p column corresponding to period ID (numeric) in \code{df}
#' @param a_pi column corresponding to child's age at beginning of time period in \code{df}
#' @param l_p column corresponding to length of time period in \code{df}
#' @param I_i column corresponding to indicator for interval censoring in \code{df}. Should be 1 if
#' child is interval-censored, 0 if right-censored.
#' @param A_i column corresponding to indicator for if a child is interval-censored across the
#' boundary of a time period in \code{df}
#' @param t_i column corresponding to age at right-censoring, if right-censored, in \code{df}
#' @param t_0i column corresponding to lower bound of interval, if interval-censored, in \code{df}
#' @param t_1i column corresponding to upper bound of interval, if interval-censored, in \code{df}
#' @param only_scale boolean for varying only the scale parameter across time period. Defaults to
#' \code{FALSE}. This option is only available for dist = "weibull"
#' @param dist distribution. Currently supports "weibull", "exponential"
#' @return A list containing: 
#' \itemize{
#' \item result: a dataframe of summarized results
#' order arrange(time). Each matrix will have rows arranged in order arrange(region).
#' \item optim: the output from \code{optim}
#' \item grad: the gradient evaluated at the MLE
#' \item variance: the finite population variance-covariance matrix
#' \item design: the survey design object
#' \item runtime: runtime for likelihood optimization
#' } 
#' 
#' @author Taylor Okonek
#' @export surv_synthetic

surv_synthetic <- function(df,
                           individual = "individual",
                           household = "household",
                           cluster = "cluster",
                           strata = "strata",
                           weights = "weights",
                           p = "p",
                           a_pi = "a_pi",
                           l_p = "l_p",
                           I_i = "I_i",
                           A_i = "A_i",
                           t_i = "t_i",
                           t_0i = "t_0i",
                           t_1i = "t_1i",
                           only_scale = FALSE,
                           dist = "weibull") {
  
  # error checking
  if (!(dist %in% c("weibull", "exponential"))) {
    stop("surv_synthetic currently only supports weibull and exponential distribution")
    if (only_scale & (dist != "weibull")) {
      stop("only_scale = TRUE is only available for dist = 'weibull'")
    }
  }
  
  # set distribution to integer
  if (dist == "weibull") {
    dist <- 1
  } else {
    dist <- 0
  }
  
  # make new df with appropriate columns
  temp <- df[,c(individual, household, cluster, strata, weights, p, a_pi, l_p,
                I_i, A_i, t_i, t_0i, t_1i)]
  df <- temp
  colnames(df) <- c("individual", "household", "cluster",
                    "strata", "weights", "p", "a_pi", "l_p",
                    "I_i", "A_i", "t_i", "t_0i", "t_1i")
  
  # make I_i == 2 if exactly observed
  df <- df %>% dplyr::mutate(I_i = ifelse((t_0i == t_1i) & I_i == 1, 2, I_i))
  
  # get number of periods
  n_periods <- length(unique(df$p))
  
  # pivot wider
  df <- df %>%
    pivot_wider(id_cols = c(individual, household, cluster, strata, weights, I_i, A_i, t_i, t_0i, t_1i),
                names_from = p,
                values_from = c(a_pi, l_p))
  
  # get column name indicators for a_pi and l_p
  a_pi_cols <- paste0("a_pi_",1:n_periods)
  l_p_cols <- paste0("l_p_",1:n_periods)
  
  # fit model
  if (only_scale) {
    message("fitting model")
    
    start_time <- Sys.time()
    optim_res <- optim(par = c(rep(-1,1), rep(15, n_periods)),
                       fn = optim_fn,
                       data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                       weights = df$weights,
                       shape_par_ids = 1,
                       dist = dist,
                       method = "BFGS",
                       hessian = TRUE)
    end_time <- Sys.time()
    end_time - start_time
    
    message("computing finite population variance")
    test_scores <- rcpp_gradient_multi(x_df = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                                       log_shapes = optim_res$par[1], 
                                       log_scales = optim_res$par[2:length(optim_res$par)], 
                                       dist = dist)
    
  } else {
    
    # weibull
    if (dist == 1) {
      message("fitting model")
      start_time <- Sys.time()
      optim_res <- optim(par = c(rep(-1,n_periods), rep(15, n_periods)),
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = 1:n_periods,
                         dist = dist,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      message("computing finite population variance")
      test_scores <- rcpp_gradient_multi(x_df = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                                         log_shapes = optim_res$par[1:n_periods], 
                                         log_scales = optim_res$par[(n_periods + 1):(n_periods * 2)], 
                                         dist = dist)
      # exponential
    } else if (dist == 0) {
      message("fitting model")
      start_time <- Sys.time()
      optim_res <- optim(par = c(rep(15, n_periods)),
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = NA,
                         dist = dist,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      message("computing finite population variance")
      test_scores <- rcpp_gradient_multi(x_df = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                                         log_shapes = 1, 
                                         log_scales = optim_res$par, 
                                         dist = 1)
      test_scores <- test_scores[,-1]
    }
    
  }
  
  # get finite pop variances
  est <- optim_res$par
  test_invinf <- solve(-optim_res$hessian)
  infl_fns <- test_scores %*% test_invinf
  design <- survey::svydesign(ids=~cluster+household, 
                              strata = ~strata,
                              weights = ~weights, 
                              data = df)
  vmat <- vcov(svytotal(infl_fns, design))
  
  # organize results to return
  if (only_scale) {
    ret_df <- data.frame(period = 1:n_periods,
                         log_shape_mean = rep(est[1], n_periods),
                         log_scale_mean = est[2:length(est)],
                         log_shape_var = rep(vmat[1,1], n_periods),
                         log_scale_var = diag(vmat)[2:length(est)])
    
    # add u5mr, nmr, imr to ret_df
    ret_df$U5MR <- p_weibull(60, log_shape = ret_df$log_shape_mean, log_scale = ret_df$log_scale_mean)
    delta_vmat <- p_weibull_deltamethod(60, par = est, vmat = vmat, shape_par_ids = 1)
    ret_df$U5MR_var <- diag(delta_vmat)
    ret_df$U5MR_upper <- ret_df$U5MR + 1.96 * sqrt(ret_df$U5MR_var)
    ret_df$U5MR_lower <- ret_df$U5MR - 1.96 * sqrt(ret_df$U5MR_var)
    
  } else {
    
    # weibull
    if (dist == 1) {
      ret_df <- data.frame(period = 1:n_periods,
                           log_shape_mean = est[1:n_periods],
                           log_scale_mean = est[(n_periods + 1):(n_periods*2)],
                           log_shape_var = diag(vmat)[1:n_periods],
                           log_scale_var = diag(vmat)[(n_periods + 1):(n_periods*2)])
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- p_weibull(60, log_shape = ret_df$log_shape_mean, log_scale = ret_df$log_scale_mean)
      delta_vmat <- p_weibull_deltamethod(60, par = est, vmat = vmat, shape_par_ids = 1:n_periods)
      ret_df$U5MR_var <- diag(delta_vmat)
      ret_df$U5MR_upper <- ret_df$U5MR + 1.96 * sqrt(ret_df$U5MR_var)
      ret_df$U5MR_lower <- ret_df$U5MR - 1.96 * sqrt(ret_df$U5MR_var)
      
      # exponential
    } else if (dist == 0) {
      ret_df <- data.frame(period = 1:n_periods,
                           log_scale_mean = est,
                           log_scale_var = diag(vmat))
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- p_weibull(60, log_shape = 1, log_scale = ret_df$log_scale_mean)
      delta_vmat <- p_exponential_deltamethod(60, par = est, vmat = vmat)
      ret_df$U5MR_var <- diag(delta_vmat)
      ret_df$U5MR_upper <- ret_df$U5MR + 1.96 * sqrt(ret_df$U5MR_var)
      ret_df$U5MR_lower <- ret_df$U5MR - 1.96 * sqrt(ret_df$U5MR_var)
    }
    
  }
  
  ret_lst <- list(result = ret_df,
                  optim = optim_res,
                  grad = test_scores,
                  variance = vmat,
                  design = design,
                  runtime = end_time - start_time)
  
  return(ret_lst)
}