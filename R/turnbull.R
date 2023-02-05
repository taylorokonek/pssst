#' Function to obtain Turbull estimates for arbitrarily censored and truncated observations
#' 
#' @param df a dataframe containing the output from \code{format_turnbull}, or optionally, dataframe
#' containing the following columns
#' @param t_0i column corresponding to lower bound of interval if interval-censored, 
#' censoring time if right-censored, observed time if exactly observed, in \code{df}
#' @param t_1i column corresponding to upper bound of interval if interval-censored, 
#' \code{Inf} if right-censored, observed time if exactly observed, in \code{df}
#' @param lefttrunc column corresponding to left-truncation time in \code{df}, 0 if not
#' left-truncated
#' @param righttrunc column corresponding to right-truncation time \code{df}, \code{Inf} if not
#' right-truncated
#' @param weights column corresponding to weights in \code{df}
#' @param jackknife_var indicator for whether or not to compute jackknife variance of the estimator.
#' Will approximate the finite sample variance when sample size is small compared to population size. Defaults
#' to \code{FALSE}. If \code{jackknife_var} = \code{TRUE}, then both cluster and strata must be specified
#' @param cluster column corresponding to cluster ID from survey design. Only used if \code{jackknife_var} = \code{TRUE}. 
#' @param strata column corresponding to strata ID from survey design. Only used if \code{jackknife_var} = \code{TRUE}. 
#' @param niter number of iterations to run the algorithm
#' @param niter_jackknife number of iterations to run each resample. Defaults to \code{niter}
#' @return a list containing the following:
#' \itemize{
#' \item a dataframe of final estimates for specific times
#' \item a matrix of estimates at each iteration, with ncol = niter
#' \item number of iterations
#' }
#' 
#' @author Taylor Okonek
#' @export turnbull
turnbull <- function(df, 
                     t_0i = "t_0i",
                     t_1i = "t_1i", 
                     lefttrunc = "lefttrunc", 
                     righttrunc = "righttrunc", 
                     weights = "weights", 
                     jackknife_var = FALSE,
                     cluster = NA,
                     strata = NA,
                     niter,
                     niter_jackknife = niter) {
  
  # error checking - TO DO
  if (jackknife_var) {
    if (is.na(cluster) | is.na(strata)) {
      stop("Must specify cluster and strata if jackknife_var = TRUE")
    }
  }
  
  # make cluster and strata factors
  df$cluster <- factor(df[,cluster])
  df$strata <- factor(df[,strata])
  
  # call rcpp_turnbull
  res <- rcpp_turnbull(niter = niter, 
                       t0 = df[,t_0i], 
                       t1 = df[,t_1i], 
                       lefttrunc = df[,lefttrunc], 
                       righttrunc = df[,righttrunc], 
                       weights = df[,weights]) 
  
  # calculate jackknife variance
  if (jackknife_var) {
    message("Calculating jackknife variance...")
    
    # get data frame containing unique clusters and stratas
    c_s_df <- df[,c(cluster, strata)] %>% dplyr::distinct()
    
    # number of strata
    h_df <- df %>% 
      dplyr::group_by(strata) %>% 
      dplyr::summarise(n_h = length(unique(cluster))) %>% 
      as.data.frame()
    h <- nrow(h_df)
    
    # number of clusters in each strata
    n_h <- h_df[,2]
    
    # number of clusters
    k <- length(unique(df[,cluster]))
    
    # calculate h_k
    h_k <- (n_h - 1)/n_h
    h_k <- rep(h_k, n_h)
    
    # loop through clusters
    r_i <- matrix(NA, nrow = nrow(res[[3]]), ncol = k)
    message(paste0("Looping through ", k, " clusters..."))
    for (i in 1:k) {
      # filter dataset to exclude cluster k
      tmp <- df[df[,cluster] != unique(df[,cluster])[i],]
      
      # adjust weights of other clusters in stratum
      which_strata <- c_s_df[,strata][c_s_df[,cluster] == i]
      other_clusts <- c_s_df[,cluster][c_s_df[,strata] == which_strata]
      other_clusts <- other_clusts[other_clusts != unique(df[,cluster])[i]]
      
      tmp[tmp[,cluster] %in% other_clusts, weights] <- n_h[which_strata]/ (n_h[which_strata] - 1)
      
      # get turnbull estimator on tiny dataset
      small_samp_est <- rcpp_turnbull(niter = niter_jackknife, 
                                      t0 = tmp[,t_0i],
                                      t1 = tmp[,t_1i],
                                      lefttrunc = tmp[,lefttrunc], 
                                      righttrunc = tmp[,righttrunc],
                                      weights = tmp[,weights])
      
      # add results to r_i
      r_i[,i] <- small_samp_est[[3]][,niter_jackknife]
    }
    
    # calculate variance
    r_i_minus_r_2 <- apply(r_i, 2, function(x) {(x - res[[3]][,niter])^2})
    for (i in 1:ncol(r_i)) {
      r_i_minus_r_2[,i] <- r_i_minus_r_2[,i] * h_k[i]
    }
    
    var_r <- rowSums(r_i_minus_r_2)
  }
  
  # create final data frame
  if (jackknife_var) {
    res_df <- data.frame(t0 = res[[1]],
                         t1 = res[[2]],
                         est = 1 - cumsum(res[[3]][,niter]),
                         var_est = var_r)
    ret_lst <- list(result = res_df,
                    iterations = apply(res[[3]],2,function(x) {1 - cumsum(x)}),
                    niter = niter,
                    niter_jackknife = niter_jackknife)
  } else {
    res_df <- data.frame(t0 = res[[1]],
                         t1 = res[[2]],
                         est = 1 - cumsum(res[[3]][,niter]))
    ret_lst <- list(result = res_df,
                    iterations = apply(res[[3]],2,function(x) {1 - cumsum(x)}),
                    niter = niter)
  }
  
  return(ret_lst)
}

