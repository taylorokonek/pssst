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
#' @param period column corresponding to period identification in \code{df}. If specified, separate
#' Turnbull estimates will be produced for each value in \code{period}.
#' @param weights column corresponding to weights in \code{df}
#' @param niter number of iterations to run the algorithm
#' @return a list containing the following:
#' \itemize{
#' \item a dataframe of final estimates for specific times (and periods)
#' \item a matrix of estimates at each iteration, with ncol = niter, or a list of such matrices for each period
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
                     period = NA,
                     weights = "weights", 
                     niter) {
  
  # Jackknife variance no longer calculated because it isn't valid
  # TO DO: delete jackknife var code (for now just comment relevant things)
  jackknife_var <- FALSE
  niter_jackknife <- niter
  # cluster <- NA
  # strata <- NA
  
  # error checking - TO DO
  if (jackknife_var) {
    if (is.na(cluster) | is.na(strata)) {
      stop("Must specify cluster and strata if jackknife_var = TRUE")
    }
  }
  
  # if there are no periods...
  if (is.na(period)) {
    # make cluster and strata factors
    # df$cluster <- factor(df[,cluster])
    # df$strata <- factor(df[,strata])
    
    # call rcpp_turnbull
    res <- rcpp_turnbull(niter = niter, 
                         t0 = df[,t_0i], 
                         t1 = df[,t_1i], 
                         lefttrunc = df[,lefttrunc], 
                         righttrunc = df[,righttrunc], 
                         weights = df[,weights],
                         set_lower = NaN,
                         set_upper = NaN) 
    
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
                                        weights = tmp[,weights],
                                        set_lower = res[[1]],
                                        set_upper = res[[2]])
        
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
    
  # if period is specified...
  } else {
    
    df$period <- df[,period]
    
    # set up list for all periods
    period_lst <- vector("list",length(unique(df$period)))
    
    # get upper and lower bounds for rcpp_turnbull in each period
    res_full <- rcpp_turnbull(niter = niter, 
                         t0 = df[,t_0i], 
                         t1 = df[,t_1i], 
                         lefttrunc = df[,lefttrunc], 
                         righttrunc = df[,righttrunc], 
                         weights = df[,weights],
                         set_lower = NaN,
                         set_upper = NaN) 
    
    # loop through periods...
    for (l in 1:length(unique(df$period))) {
      # filter dataframe
      df_p <- df %>% dplyr::filter(period == sort(unique(df$period))[l])
      
      # make cluster and strata factors
      # df_p$cluster <- factor(df_p[,cluster])
      # df_p$strata <- factor(df_p[,strata])
      
      # call rcpp_turnbull
      res <- rcpp_turnbull(niter = niter, 
                           t0 = df_p[,t_0i], 
                           t1 = df_p[,t_1i], 
                           lefttrunc = df_p[,lefttrunc], 
                           righttrunc = df_p[,righttrunc], 
                           weights = df_p[,weights],
                           set_lower = res_full[[1]],
                           set_upper = res_full[[2]]) 
      
      # calculate jackknife variance
      if (jackknife_var) {
        message("Calculating jackknife variance...")
        
        # get data frame containing unique clusters and stratas
        c_s_df <- df_p[,c(cluster, strata)] %>% dplyr::distinct()
        
        # number of strata
        h_df <- df_p %>% 
          dplyr::group_by(strata) %>% 
          dplyr::summarise(n_h = length(unique(cluster))) %>% 
          as.data.frame()
        h <- nrow(h_df)
        
        # number of clusters in each strata
        n_h <- h_df[,2]
        
        # number of clusters
        k <- length(unique(df_p[,cluster]))
        
        # calculate h_k
        h_k <- (n_h - 1)/n_h
        h_k <- rep(h_k, n_h)
        
        # loop through clusters
        r_i <- matrix(NA, nrow = nrow(res[[3]]), ncol = k)
        message(paste0("Looping through ", k, " clusters in period ", sort(unique(df$period))[l], "..."))
        
        for (i in 1:k) {
          
          which_cluster <- sort(unique(df_p[,cluster]))[i]
          # filter dataset to exclude cluster k
          tmp <- df_p[df_p[,cluster] != which_cluster,]
          
          # adjust weights of other clusters in stratum
          which_strata <- c_s_df[,strata][c_s_df[,cluster] == which_cluster]
          other_clusts <- c_s_df[,cluster][c_s_df[,strata] == which_strata]
          other_clusts <- other_clusts[other_clusts != which_cluster]
          
          # if only one PSU in strata, cannot calculate variance
          if (n_h[which_strata] == 1) {
            message(paste0("Stratum ", which_strata, " has only one PSU"))
            next
          }
          
          tmp[tmp[,cluster] %in% other_clusts, weights] <- n_h[which_strata]/ (n_h[which_strata] - 1)
          
          # get turnbull estimator on tiny dataset
          small_samp_est <- rcpp_turnbull(niter = niter_jackknife, 
                                          t0 = tmp[,t_0i],
                                          t1 = tmp[,t_1i],
                                          lefttrunc = tmp[,lefttrunc], 
                                          righttrunc = tmp[,righttrunc],
                                          weights = tmp[,weights],
                                          set_lower = res_full[[1]],
                                          set_upper = res_full[[2]])
          
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
      
      # add results to existing (outside of period loop)
      period_lst[[l]] <- ret_lst
    }
    
    # combine results across periods
    if (jackknife_var) {
      
      # combine results from all periods into one long dataframe
      temp_lst <- vector("list", length = length(period_lst))
      temp_lst_iter <- vector("list", length = length(period_lst))
      for (i in 1:length(period_lst)) {
        period_lst[[i]][[1]][,period] <- sort(unique(df$period))[i]
        temp_lst[[i]] <- period_lst[[i]][[1]][,c(period, "t0", "t1", "est", "var_est")]
        temp_lst_iter[[i]] <- period_lst[[i]][[2]]
      }
      res_df <- do.call(rbind, temp_lst)
      
      ret_lst <- list(result = res_df,
                      iterations = temp_lst_iter,
                      niter = niter,
                      niter_jackknife = niter_jackknife)
    } else {
      
      # combine results from all periods into one long dataframe
      temp_lst <- vector("list", length = length(period_lst))
      temp_lst_iter <- vector("list", length = length(period_lst))
      for (i in 1:length(period_lst)) {
        period_lst[[i]][[1]][,period] <- sort(unique(df$period))[i]
        temp_lst[[i]] <- period_lst[[i]][[1]][,c(period, "t0", "t1", "est")]
        temp_lst_iter[[i]] <- period_lst[[i]][[2]]
      }
      res_df <- do.call(rbind, temp_lst)
      
      ret_lst <- list(result = res_df,
                      iterations = temp_lst_iter,
                      niter = niter)
    }
    
  }
  
  return(ret_lst)
}

