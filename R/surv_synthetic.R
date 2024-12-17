#' Fit a parametric survival model to DHS data across multiple time periods.
#' 
#' Fits a parametric survival model to DHS data (formatted from \code{format_dhs}) to get survey-weighted
#' parameter estimates and finite population variances for parametric survival curves in each time period.
#' 
#' We use the same parameterization of the generalized gamma distribution
#' as used in \code{\link{flexsurv}}, as opposed to the original parameterization 
#' described in Stacy 1962, where if \eqn{\omega \sim Gamma(k, 1)}, then \eqn{x = \exp(\omega / shape + \log(scale))}
#' follows the original generalized gamma distribution. With \eqn{shape = b > 0} and \eqn{scale = a > 0}, \eqn{x}
#' has the pdf
#' 
#' \deqn{f(x \mid a, b, k) = \frac{b}{\Gamma(k)} \frac{x^{bk - 1}}{a^{bk}} \exp(-(x/a)^b)}
#' 
#' The alternative parameterization developed in Prentice (1974) is typically preferred,
#' as it is more numerically stable in some cases. Under this alternative, if 
#' \eqn{\gamma \sim Gamma(Q^{-2}, 1)} and \eqn{\omega = \log(Q^2 \gamma) / Q}, then \eqn{x = \exp(\mu + \sigma \omega)}
#' follows the generalized gamma distribution with pdf
#' 
#' \deqn{f(x \mid \mu, \sigma, Q) = \frac{|Q| (Q^{-2})^{Q^{-2}}}{\sigma x \Gamma(Q^{-2})} \exp(Q^{-2}(Q \omega - \exp(Q \omega)))}
#' 
#' The "etsp" distribution uses the exponentially-trunacted shifted power family of hazards defined by
#' 
#' \deqn{h(x) = a(x + c)^{-p} e^{-b x}}
#' 
#' in Scholey (2019), but sets \eqn{c = 0}.
#' 
#' @param df a dataframe containing the output from \code{format_dhs}, or optionally, dataframe
#' containing the following columns
#' @param individual column corresponding to individual ID in \code{df}
#' @param survey boolean for whether or not your data comes from a survey with a specified probability design.
#' If true, must specify household, cluster, strata, and weights parameters in the function, and a finite
#' population variance based on a pseudo-likelihood will be returned for your survey-weighted estimates. 
#' If false, superpopulation estimates will be returned with the usual asymptotic variance estimator.
#' Defaults to \code{TRUE}.
#' @param household column corresponding to household ID in \code{df}
#' @param cluster column corresponding to cluster ID in \code{df}
#' @param strata column corresponding to strata ID in \code{df}
#' @param weights column corresponding to weights in \code{df}. If survey = FALSE, you may
#' still specify a weights column, but note that the superpopulation variance will be returned
#' as opposed to the finite population variance. 
#' @param p column corresponding to period ID (must be integer-valued, numeric) in \code{df}
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
#' @param numerical_grad boolean for whether gradient should be calculated numerically or
#' analytically. Analytical gradient is faster, but only available for Weibull and Exponential distributions
#' at the moment.
#' @param dist distribution. Currently supports "weibull", "exponential", 
#' "piecewise_exponential", "gengamma", "lognormal", "gompertz", "etsp" (exponentially-truncated shifted power family),
#' "loglogistic", "dagum"
#' @param breakpoints if distribution is "piecewise_exponential", the breakpoints (in months) where
#' the distribution should be divided
#' @param init_vals an optional vector of initial values at which to start the optimizer for the 
#' parameters. Must specify the appropriate number of parameters for the given distribution /
#' number of periods.
#' @param etsp_c If dist is etsp, fixed value to use for the c parameter. Default is zero.
#' @return A list containing: 
#' \itemize{
#' \item result: a dataframe of summarized results
#' order arrange(time). Each matrix will have rows arranged in order arrange(region).
#' \item optim: the output from \code{optim}
#' \item grad: the gradient evaluated at the MLE
#' \item variance: the finite population variance-covariance matrix
#' \item design: the survey design object
#' \item runtime: runtime for likelihood optimization
#' \item initial_values: initial values for parameters used in likelihood maximization
#' } 
#' 
#' @author Taylor Okonek
#' 
#' @references Stacy, E. W. (1962). A generalization of the gamma
##' distribution.  Annals of Mathematical Statistics 33:1187-92.
##' 
##' Prentice, R. L. (1974). A log gamma model and its maximum likelihood
##' estimation. Biometrika 61(3):539-544.
##' 
##' Scholey, J (2019). The Age-Trajectory of Infant Mortality in the United States: 
##' Parametric Models and Generative Mechanisms. Annual meeting of the Population 
##' Association of America, Austin, TX. 2019.
##' 
#' @export surv_synthetic

surv_synthetic <- function(df,
                           individual = "individual",
                           survey = TRUE,
                           household = NULL,
                           cluster = NULL,
                           strata = NULL,
                           weights = NULL,
                           p = "p",
                           a_pi = "a_pi",
                           l_p = "l_p",
                           I_i = "I_i",
                           A_i = "A_i",
                           t_i = "t_i",
                           t_0i = "t_0i",
                           t_1i = "t_1i",
                           only_scale = FALSE,
                           numerical_grad = FALSE,
                           dist = "weibull",
                           breakpoints = NA,
                           init_vals = NA,
                           etsp_c = 0) {
  
  # error checking
  if (survey & is.null(weights)) {
    stop("If survey == true, must specify at least weights columns")
  }
  
  # format ids = ~household+cluster
  if (!is.null(household)) {
    if (!is.null(cluster)) {
      ids_form <- formula(paste0("~",cluster,"+",household))
    } else {
      ids_form <- formula(paste0("~",household))
    }
  } else {
    if (!is.null(cluster)) {
      ids_form <- formula(paste0("~",cluster))
    } else {
      ids_form <- formula("~1")
    }
  }
  
  # format strata = ~strata
  if (!is.null(strata)) {
    strata_form <- formula(paste0("~",strata))
  } else {
    strata_form = NULL
  }
  
  # format weights = ~weights
  if (!is.null(weights)) {
    weights_form <- formula(paste0("~",weights))
  } else {
    weights_form <- NULL
  }
  
  # distribution errors
  if (!(dist %in% c("weibull", "exponential", "piecewise_exponential", 
                    "gengamma", "lognormal", "gompertz","etsp", "loglogistic", 
                    "dagum"))) {
    stop("distribution not currently supported in surv_synthetic")
  }
  if (only_scale & (dist != "weibull")) {
    stop("only_scale = TRUE is only available for dist = 'weibull'")
  }
  if ((dist == "piecewise_exponential") & is.na(breakpoints[1])) {
    stop("breakpoints must be set if using piecewise exponential distribution")
  }
  if (!is.na(breakpoints[1]) & dist != "piecewise_exponential") {
    message(paste0("breakpoints not available for distribution ",dist,", and will not be used"))
  }
  
  if (!(dist %in% c("exponential", "weibull"))) {
    if (!numerical_grad) {
      stop("analytical gradient only available for exponential and weibull dist. must set numerical_grad = TRUE")
    }
  }
  
  # Check that t_0i <= t_1i everywhere
  if(sum(df[, t_0i] > df[, t_1i], na.rm=T) != 0) {
    stop("at least one observation has an interval-censored time where t_0 > t_1")
  }
  
  # remove 0 and Inf from breakpoints, if needed, and sort
  if (!is.na(breakpoints[1])) {
    breakpoints <- sort(breakpoints[!(breakpoints %in% c(0, Inf))])
  }
  
  # set distribution to integer
  if (dist == "weibull") {
    dist <- 1
  } else if (dist == "exponential") {
    dist <- 0
  } else if (dist == "piecewise_exponential") {
    dist <- 2
  } else if (dist == "gengamma") {
    dist <- 3
  } else if (dist == "lognormal") {
    dist <- 4 
  } else if (dist == "gompertz") {
    dist <- 5 
  } else if (dist == "etsp") {
    dist <- 6 
  } else if (dist == "loglogistic") {
    dist <- 7 # loglogistic
  } else if (dist == "dagum") {
    dist <- 8 
  }
  # make new df with appropriate columns
  temp <- df[,c(individual, household, cluster, strata, weights, p, a_pi, l_p,
                I_i, A_i, t_i, t_0i, t_1i)]
  df <- temp
  survey_cols_include <- !c(is.null(household), is.null(cluster), is.null(strata), is.null(weights))
  colnames(df) <- c("individual", c("household", "cluster","strata","weights")[which(survey_cols_include)],
                    "p", "a_pi", "l_p",
                    "I_i", "A_i", "t_i", "t_0i", "t_1i")
  
  # make I_i == 2 if exactly observed
  df <- df %>% dplyr::mutate(I_i = ifelse((t_0i == t_1i) & I_i == 1, 2, I_i))
  
  # get number of periods
  n_periods <- length(unique(df$p))
  period_names <- unique(df$p)
  
  # Throw error if time period is neither an integer nor a string
  if (!(class(df$p) %in% c("numeric","integer"))) {
    stop("Period must be an integer-valued, numeric column.")
  }
  if (sum(round(df$p) - df$p) != 0) {
    stop("Period must be an integer-valued, numeric column.")
  }
  
  # Throw warning if more than 20 time periods are present
  if (n_periods > 20) {
    warning(paste0("This warning appears if you have more than 20 time periods present in your data. ", n_periods, " periods are present in your data. Please make sure you specified the correct column for time period within the function."))
  }
  
  # check that a_pi is unique for each individual in each time period
  num_diff_a_pis <- df %>% group_by(individual) %>%
    summarize(temp = length(unique(a_pi))) %>%
    dplyr::select(temp) %>% unlist %>% unname() %>% unique() %>% length()
  if (num_diff_a_pis > 1) {
    stop("Each individual much have a distinct a_pi value in each time period.")
  }
  
  # check that init_vals is the correct length (if specified)
  if (!is.na(init_vals[1])) {
    if (dist %in% c(4,5,7)) { # lognormal, gompertz
      if (length(init_vals) != (2 * n_periods)) {
        stop("Incorrect number of initial values supplied")
      }
    } else if (dist == 0) { # exponential
      if (length(init_vals) != (n_periods)) {
        stop("Incorrect number of initial values supplied")
      }
    } else if (dist == 1) { # weibull
      if (only_scale) {
        if (length(init_vals) != (1 + n_periods)) {
          stop("Incorrect number of initial values supplied")
        }
      } else {
        if (length(init_vals) != (2 * n_periods)) {
          stop("Incorrect number of initial values supplied")
        }
      }
    } else if (dist %in% c(3,6)) {# gengamma, etsp
      if (length(init_vals) != (3 * n_periods)) {
        stop("Incorrect number of initial values supplied")
      }
    } else { # piecewise exponential (dist == 2)
      if (length(init_vals) != (n_periods * (length(breakpoints) + 1))) {
        stop("Incorrect number of initial values supplied")
      }
    }
  }
  
  # pivot wider
  df <- df %>%
    tidyr::pivot_wider(id_cols = c(individual, household, cluster, strata, weights, I_i, A_i, t_i, t_0i, t_1i),
                       names_from = p,
                       values_from = c(a_pi, l_p))
  
  # get column name indicators for a_pi and l_p
  a_pi_cols <- paste0("a_pi_",period_names)
  l_p_cols <- paste0("l_p_",period_names)

 
  
  # if no survey design and weights column unspecified, set weights = 1 for all individuals
  if (!survey & is.null(weights)) {
    df$weights <- 1
  }
  
  # fit model
  if (only_scale) {
    message("fitting model")
    
    if (is.na(init_vals[1])) {
      init_vals <- c(rep(-1,1), rep(15, n_periods))
    }
    
    start_time <- Sys.time()
    optim_res <- optim(par = init_vals,
                       fn = optim_fn,
                       data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                       weights = df$weights,
                       shape_par_ids = 1,
                       dist = dist,
                       breakpoints = breakpoints,
                       num_periods = n_periods,
                       method = "BFGS",
                       hessian = TRUE)
    end_time <- Sys.time()
    end_time - start_time
    
    if (survey) {
      message("computing finite population variance")
      if (numerical_grad) {
        test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
        for (i in 1:nrow(df)) {
          test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                            x = optim_res$par, 
                                            data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                            weights = 1,
                                            shape_par_ids = 1,
                                            dist = dist,
                                            breakpoints = breakpoints,
                                            num_periods = n_periods,
                                            etsp_c = etsp_c)
        }
      } else {
        test_scores <- rcpp_gradient_multi(x_df = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                                           log_shapes = optim_res$par[1], 
                                           log_scales = optim_res$par[2:length(optim_res$par)], 
                                           dist = dist)
      }
    }
    
  } else {
    
    # weibull
    if (dist == 1) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(-1,n_periods), rep(15, n_periods))
      }

      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = 1:n_periods,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        if (numerical_grad) {
          test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
          for (i in 1:nrow(df)) {
            test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                              x = optim_res$par, 
                                              data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                              weights = 1,
                                              shape_par_ids = 1:n_periods,
                                              dist = dist,
                                              breakpoints = breakpoints,
                                              num_periods = n_periods)
          }
        } else {
          test_scores <- rcpp_gradient_multi(x_df = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                                             log_shapes = optim_res$par[1:n_periods], 
                                             log_scales = optim_res$par[(n_periods + 1):(n_periods * 2)], 
                                             dist = dist)
        }
      }
      
      
      # exponential
    } else if (dist == 0) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(15, n_periods))
      }

      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = NA,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        if (numerical_grad) {
          test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
          for (i in 1:nrow(df)) {
            test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                              x = optim_res$par, 
                                              data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                              weights = 1,
                                              shape_par_ids = NA,
                                              dist = 1,
                                              breakpoints = breakpoints,
                                              num_periods = n_periods)
            test_scores <- test_scores[,-1]
          }
        } else {
          test_scores <- rcpp_gradient_multi(x_df = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                                             log_shapes = 1, 
                                             log_scales = optim_res$par, 
                                             dist = 1)
          test_scores <- test_scores[,-1]
        }
      }
      
      # piecewise exponential
    } else if (dist == 2) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(15, n_periods * (length(breakpoints) + 1)))
      }
      
      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = NA,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        
        test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
        for (i in 1:nrow(df)) {
          test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                            x = optim_res$par, 
                                            data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                            weights = 1,
                                            shape_par_ids = NA,
                                            dist = dist,
                                            breakpoints = breakpoints,
                                            num_periods = n_periods)
        }
      }
      
      # generalized gamma
    } else if (dist == 3) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(1,n_periods), rep(1, n_periods * 2))
      }
      
      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = 1:n_periods,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        
        test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
        for (i in 1:nrow(df)) {
          test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                            x = optim_res$par, 
                                            data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                            weights = 1,
                                            shape_par_ids = 1:n_periods,
                                            dist = dist,
                                            breakpoints = breakpoints,
                                            num_periods = n_periods)
        }
      }
      
      # lognormal
    } else if (dist == 4) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(1.1, n_periods), rep(-2, n_periods))
      }
      
      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = 1:n_periods,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE)
      # glimpse(optim_res)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        
        test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
        for (i in 1:nrow(df)) {
          test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                            x = optim_res$par, 
                                            data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                            weights = 1,
                                            shape_par_ids = 1:n_periods,
                                            dist = dist,
                                            breakpoints = breakpoints,
                                            num_periods = n_periods)
        }
      }
      
      # gompertz
    } else if (dist == 5) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(-10, n_periods ), rep(4, n_periods))
      }
      
      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = 1:n_periods,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        
        test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
        for (i in 1:nrow(df)) {
          test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                            x = optim_res$par, 
                                            data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                            weights = 1,
                                            shape_par_ids = 1:n_periods,
                                            dist = dist,
                                            breakpoints = breakpoints,
                                            num_periods = n_periods)
        }
      } 
      
      # etsp
    } else if (dist == 6) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(-4, n_periods), rep(-4, n_periods * 2))
      }
      
      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = 1:n_periods,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE,
                         etsp_c = etsp_c)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        
        test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
        for (i in 1:nrow(df)) {
          test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                            x = optim_res$par, 
                                            data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                            weights = 1,
                                            shape_par_ids = 1:n_periods,
                                            dist = dist,
                                            breakpoints = breakpoints,
                                            num_periods = n_periods,
                                            etsp_c = etsp_c)
        }
      }
      
      # loglogistc
    } else if (dist == 7) {
      message("fitting model")
      
      if (is.na(init_vals[1])) {
        init_vals <- c(rep(-4, n_periods), rep(-4, n_periods))
      }
      
      start_time <- Sys.time()
      optim_res <- optim(par = init_vals,
                         fn = optim_fn,
                         data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                         weights = df$weights,
                         shape_par_ids = 1:n_periods,
                         dist = dist,
                         breakpoints = breakpoints,
                         num_periods = n_periods,
                         method = "BFGS",
                         hessian = TRUE)
      end_time <- Sys.time()
      end_time - start_time
      
      if (survey) {
        message("computing finite population variance")
        
        test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
        for (i in 1:nrow(df)) {
          test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                            x = optim_res$par, 
                                            data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                            weights = 1,
                                            shape_par_ids = 1:n_periods,
                                            dist = dist,
                                            breakpoints = breakpoints,
                                            num_periods = n_periods)
        }
      }
    } else if (dist == 8) {
    message("fitting model")
    
    if (is.na(init_vals[1])) {
      init_vals <- c(rep(-4, n_periods), rep(-4, n_periods), rep(-4, n_periods))
    }

    start_time <- Sys.time()
    optim_res <- optim(par = init_vals,
                       fn = optim_fn,
                       data = df[,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)] %>% as.data.frame(),
                       weights = df$weights,
                       shape_par_ids = 1:n_periods,
                       # shape_par_ids = 1:n_periods,
                       dist = dist,
                       breakpoints = breakpoints,
                       num_periods = n_periods,
                       method = "BFGS",
                       hessian = TRUE)
    end_time <- Sys.time()
    end_time - start_time
    
    if (survey) {
      message("computing finite population variance")
      
      test_scores <- matrix(nrow = nrow(df), ncol = length(optim_res$par))
      for (i in 1:nrow(df)) {
        test_scores[i,] <- numDeriv::grad(optim_fn_grad, 
                                          x = optim_res$par, 
                                          data = df[i,c("I_i","A_i","t_i","t_0i","t_1i", a_pi_cols, l_p_cols)],
                                          weights = 1,
                                          shape_par_ids = 1:n_periods,
                                          dist = dist,
                                          breakpoints = breakpoints,
                                          num_periods = n_periods)
      }
    }
  }
  }
  
  # get finite pop variances
  est <- optim_res$par
  
  if (survey) {
    test_invinf <- solve(-optim_res$hessian)
    infl_fns <- test_scores %*% test_invinf
    design <- survey::svydesign(ids = ids_form, 
                                strata = strata_form,
                                weights = weights_form, 
                                data = df)
    vmat <- vcov(survey::svytotal(infl_fns, design))
  } else {
    message("computing superpopulation variance")
    vmat <- solve(optim_res$hessian)
  }
  
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
    
    ret_df$IMR <- p_weibull(12, log_shape = ret_df$log_shape_mean, log_scale = ret_df$log_scale_mean)
    ret_df$NMR <- p_weibull(1, log_shape = ret_df$log_shape_mean, log_scale = ret_df$log_scale_mean)
    
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
      
      ret_df$IMR <- p_weibull(12, log_shape = ret_df$log_shape_mean, log_scale = ret_df$log_scale_mean)
      ret_df$NMR <- p_weibull(1, log_shape = ret_df$log_shape_mean, log_scale = ret_df$log_scale_mean)
      
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
      
      ret_df$IMR <- p_weibull(12, log_shape = 1, log_scale = ret_df$log_scale_mean)
      ret_df$NMR <- p_weibull(1, log_shape = 1, log_scale = ret_df$log_scale_mean)
      
      # piecewise exponential
    } else if (dist == 2) {
      # get breakpoints in form (a,b]
      breakpoints_full <- sort(unique(c(0,breakpoints,Inf)))
      breakpoints_lower <- breakpoints_full[-length(breakpoints_full)]
      breakpoints_upper <- breakpoints_full[-1]
      breakpoint_intervals <- paste0("[", breakpoints_lower, ",", breakpoints_upper, "]")
      
      ret_df <- data.frame(period = 1:n_periods)
      ret_df <- cbind(ret_df, data.frame(matrix(est, nrow = n_periods, byrow = TRUE)))
      colnames(ret_df)[2:ncol(ret_df)] <- paste0("log_scale_mean_",breakpoint_intervals)
      
      ret_df <- cbind(ret_df, data.frame(matrix(diag(vmat), nrow = n_periods, byrow = TRUE)))
      colnames(ret_df)[(3 + length(breakpoints)):ncol(ret_df)] <- paste0("log_scale_var_", breakpoint_intervals)
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- NA
      ret_df$IMR <- NA
      ret_df$NMR <- NA
      for (i in 1:nrow(ret_df)) {
        ret_df$U5MR[i] <- p_piecewise_exponential(60, log_scales = unlist(ret_df[i,c(2:(length(breakpoints) + 2))]), breakpoints = breakpoints)
        ret_df$IMR[i] <- p_piecewise_exponential(12, log_scales = unlist(ret_df[i,c(2:(length(breakpoints) + 2))]), breakpoints = breakpoints)
        ret_df$NMR[i] <- p_piecewise_exponential(1, log_scales = unlist(ret_df[i,c(2:(length(breakpoints) + 2))]), breakpoints = breakpoints)
      }
      
      # generalized gamma
    } else if (dist == 3) {
      ret_df <- data.frame(period = 1:n_periods,
                           log_sigma_mean = est[1:n_periods])
      ret_df <- cbind(ret_df, data.frame(est[(n_periods +1 ):length(est)] %>% matrix(nrow = n_periods, byrow = TRUE)))
      colnames(ret_df)[c(3:4)] <- c("log_mu_mean", "log_Q_mean")
      
      ret_df$log_sigma_var <- diag(vmat)[1:n_periods]
      ret_df <- cbind(ret_df, data.frame(diag(vmat)[(n_periods +1 ):length(est)] %>% matrix(nrow = n_periods, byrow = TRUE)))
      colnames(ret_df)[c(6:7)] <- c("log_mu_var", "log_Q_var")
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- NA
      ret_df$IMR <- NA
      ret_df$NMR <- NA
      for (i in 1:nrow(ret_df)) {
        ret_df$U5MR[i] <- rcpp_F_gengamma(60, exp(ret_df$log_Q_mean)[i], exp(ret_df$log_mu_mean)[i], exp(ret_df$log_sigma_mean)[i], 1, 0)
        ret_df$IMR[i] <- rcpp_F_gengamma(12, exp(ret_df$log_Q_mean)[i], exp(ret_df$log_mu_mean)[i], exp(ret_df$log_sigma_mean)[i], 1, 0)
        ret_df$NMR[i] <- rcpp_F_gengamma(1, exp(ret_df$log_Q_mean)[i], exp(ret_df$log_mu_mean)[i], exp(ret_df$log_sigma_mean)[i], 1, 0)
      }
      
      # lognormal
    } else if (dist == 4) {
      ret_df <- data.frame(period = 1:n_periods,
                           log_sigma_mean = est[1:n_periods],
                           log_1overmu_mean = est[(n_periods + 1):(n_periods*2)],
                           mu_mean = 1/exp(est[(n_periods + 1):(n_periods*2)]),
                           log_sigma_var = diag(vmat)[1:n_periods],
                           log_1overmu_var = diag(vmat)[(n_periods + 1):(n_periods*2)])
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- plnorm(60, meanlog = ret_df$mu_mean, sdlog = exp(ret_df$log_sigma_mean))
      ret_df$IMR <- plnorm(12, meanlog = ret_df$mu_mean, sdlog = exp(ret_df$log_sigma_mean))
      ret_df$NMR <- plnorm(1, meanlog = ret_df$mu_mean, sdlog = exp(ret_df$log_sigma_mean))
      
      # gompertz
    } else if (dist == 5) {
      ret_df <- data.frame(period = 1:n_periods,
                           log_shape_mean = est[1:n_periods],
                           log_scale_mean = est[(n_periods + 1):(n_periods*2)],
                           log_shape_var = diag(vmat)[1:n_periods],
                           log_scale_var = diag(vmat)[(n_periods + 1):(n_periods*2)])
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- NA
      ret_df$IMR <- NA
      ret_df$NMR <- NA
      for (i in 1:nrow(ret_df)) {
        ret_df$U5MR[i] <- rcpp_F_gompertz(60, rate = 1/exp(ret_df$log_scale_mean[i]), shape = exp(ret_df$log_shape_mean[i]), 1, 0)
        ret_df$IMR[i] <- rcpp_F_gompertz(12, rate = 1/exp(ret_df$log_scale_mean[i]), shape = exp(ret_df$log_shape_mean[i]), 1, 0)
        ret_df$NMR[i] <- rcpp_F_gompertz(1, rate = 1/exp(ret_df$log_scale_mean[i]), shape = exp(ret_df$log_shape_mean[i]), 1, 0)
      }
      
      # etsp
    } else if (dist == 6) {
      temp_mat <- est[(n_periods + 1):length(est)] %>% matrix(ncol = 2, byrow = TRUE)
      temp_vmat <- diag(vmat)[(n_periods + 1):length(est)] %>% matrix(ncol = 2, byrow = TRUE)
      
      ret_df <- data.frame(period = 1:n_periods,
                           log_a_mean = est[1:n_periods],
                           log_b_mean = temp_mat[,1],
                           log_p_mean = temp_mat[,2],
                           log_a_var = diag(vmat)[1:n_periods],
                           log_b_var = temp_vmat[,1],
                           log_p_var = temp_vmat[,2])
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- NA
      ret_df$IMR <- NA
      ret_df$NMR <- NA
      for (i in 1:nrow(ret_df)) {
        MR <- 1 - exp(-H_etsp(0, c(1, 12, 60),
                               a = exp(ret_df$log_a_mean[i]), 
                               b = exp(ret_df$log_b_mean[i]), 
                               c = etsp_c,
                               p = exp(ret_df$log_p_mean[i])))
        ret_df$NMR[i] <- MR[1]
        ret_df$IMR[i] <- MR[2]
        ret_df$IMR[i] <- MR[3]
      }
      
    } else if (dist == 7) {
      ret_df <- data.frame(period = 1:n_periods,
                           log_shape_mean = est[1:n_periods],
                           log_scale_mean = est[(n_periods + 1):(n_periods*2)],
                           log_shape_var = diag(vmat)[1:n_periods],
                           log_scale_var = diag(vmat)[(n_periods + 1):(n_periods*2)])
      
      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- NA
      ret_df$NMR <- NA
      ret_df$IMR <- NA
      for (i in 1:nrow(ret_df)) {
        ret_df$U5MR[i] <- rcpp_F_loglogistic(60, shape = exp(ret_df$log_shape_mean[i]), scale = exp(ret_df$log_scale_mean[i]), 1, 0)
        ret_df$NMR[i] <- rcpp_F_loglogistic(1, shape = exp(ret_df$log_shape_mean[i]), scale = exp(ret_df$log_scale_mean[i]), 1, 0)
        ret_df$IMR[i] <- rcpp_F_loglogistic(12, shape = exp(ret_df$log_shape_mean[i]), scale = exp(ret_df$log_scale_mean[i]), 1, 0)
      }
    }
    else if (dist == 8) {

      ret_df <- data.frame(period = 1:n_periods,
                           log_p_mean = est[1:n_periods])
      
      ret_df <- cbind(ret_df, data.frame(est[(n_periods+1):length(est)] %>%
                                           matrix(nrow = n_periods, byrow = TRUE)))
      
      colnames(ret_df)[c(3:4)] <- c("log_a_mean", "log_b_mean")
      
      ret_df$log_p_var <- diag(vmat)[1:n_periods]
      ret_df <- cbind(ret_df, data.frame(diag(vmat)[(n_periods+1):length(est)] %>%
                                           matrix(nrow = n_periods, byrow = TRUE)))
      colnames(ret_df)[c(6:7)] <- c("log_a_var", "log_b_var")
      
      
      

      # add u5mr, nmr, imr to ret_df
      ret_df$U5MR <- NA
      ret_df$NMR <- NA
      ret_df$IMR <- NA
      
      for (i in 1:nrow(ret_df)) {
        ret_df$U5MR[i] <- rcpp_F_dagum(60,
                                       shape1 = exp(ret_df$log_p_mean[i]),
                                       shape2 = exp(ret_df$log_a_mean[i]), 
                                       scale = exp(ret_df$log_b_mean[i]), 1, 0)
        ret_df$NMR[i] <- rcpp_F_dagum(1,
                                      shape1 = exp(ret_df$log_p_mean[i]),  
                                      shape2 = exp(ret_df$log_a_mean[i]), 
                                      scale = exp(ret_df$log_b_mean[i]), 1, 0)
        ret_df$IMR[i] <- rcpp_F_dagum(12, 
                                      shape1 = exp(ret_df$log_p_mean[i]),  
                                      shape2 = exp(ret_df$log_a_mean[i]), 
                                      scale = exp(ret_df$log_b_mean[i]), 1, 0)
      }
    
    }
  }
  
  if (survey) {
    ret_lst <- list(result = ret_df,
                    optim = optim_res,
                    grad = test_scores,
                    variance = vmat,
                    design = design,
                    runtime = end_time - start_time,
                    initial_values = init_vals)
  } else {
    ret_lst <- list(result = ret_df,
                    optim = optim_res,
                    variance = vmat,
                    runtime = end_time - start_time,
                    initial_values = init_vals)
  }
  
  return(ret_lst)
}