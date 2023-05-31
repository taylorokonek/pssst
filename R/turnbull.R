#' Function to obtain Turbull estimates for arbitrarily censored and truncated 
#' observations
#' 
#' @param df a dataframe containing the output from \code{format_turnbull}, or 
#' optionally, dataframe containing the following columns
#' @param t_0i column corresponding to lower bound of interval if 
#' interval-censored, censoring time if right-censored, observed time if 
#' exactly observed, in \code{df}
#' @param t_1i column corresponding to upper bound of interval if 
#' interval-censored, \code{Inf} if right-censored, observed time if exactly 
#' observed, in \code{df}
#' @param lefttrunc column corresponding to left-truncation time in \code{df}, 
#' 0 if not left-truncated
#' @param righttrunc column corresponding to right-truncation time \code{df}, 
#' \code{Inf} if not right-truncated
#' @param period column corresponding to period identification in \code{df}. 
#' If specified, separate Turnbull estimates will be produced for each value 
#' in \code{period}.
#' @param weights column corresponding to weights in \code{df}
#' @param bootstrap_var boolean indicating whether or not to calculate a finite
#' population bootstrap variance for the Turnbull estimator. If 
#' \code{bootstrap_var = TRUE}, both \code{cluster} and \code{strata} must be 
#' specified.
#' @param niter number of iterations to run the algorithm
#' @param niter_bootstrap number of bootstrap samples to take when calculating
#' bootstrap variance. If not specified, defaults to \code{niter}.
#' @return a list containing the following:
#' \itemize{
#' \item a dataframe of final estimates for specific times (and periods)
#' \item a matrix of estimates at each iteration, with ncol = niter, or a list 
#' of such matrices for each period
#' \item if \code{bootstrap_var = TRUE}, a list containing a matrix of 
#' bootstrap samples for each period. If only one time period, just a single
#' matrix
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
                     bootstrap_var = FALSE,
                     cluster = NULL,
                     strata = NULL,
                     niter,
                     niter_bootstrap = niter) {
  
  # error checking - TO DO
  if (bootstrap_var) {
    if (is.null(cluster) | is.null(strata)) {
      stop("Must specify cluster and strata if bootstrap_var = TRUE")
    }
  }
  
  # if there are no periods...
  if (is.na(period)) {
    # make cluster and strata factors
    if (!is.null(cluster)) {
      df$cluster <- factor(df[,cluster])
    }
    if (!is.null(strata)) {
      df$strata <- factor(df[,strata])
    }
    
    # call rcpp_turnbull
    res <- rcpp_turnbull(niter = niter, 
                         t0 = df[,t_0i], 
                         t1 = df[,t_1i], 
                         lefttrunc = df[,lefttrunc], 
                         righttrunc = df[,righttrunc], 
                         weights = df[,weights],
                         set_lower = NaN,
                         set_upper = NaN) 
    
    # calculate bootstrap variance
    if (bootstrap_var) {
      message("Calculating bootstrap variance...")
      
      # get data frame containing unique clusters and stratas
      c_s_df <- df[,c(cluster, strata)] %>% dplyr::distinct()
      
      # number of strata
      h_df <- df %>% 
        dplyr::group_by(strata) %>% 
        dplyr::summarise(n_h = length(unique(cluster))) %>% 
        as.data.frame()
      h <- nrow(h_df)
      
      # vector of unique clusters
      unique_clusts <- unique(df[,cluster])
      
      # set up matrix for bootstrap results
      boot_mat <- matrix(NA, nrow = length(res_full[[1]]), ncol = niter_bootstrap)
      
      # for i in 1:niter_bootstrap...
      # - sample clusters with replacement
      # - multiply design weights by correction factor n_h / (n_h - 1)
      # - calculate turnbull estimator
      
      for (i in 1:niter_bootstrap) {
        # sample clusters with replacement within each strata
        clust_samp <- c()
        for (i in 1:length(unique(c_s_df$strata))) {
          temp_clusts <- c_s_df[c_s_df[,strata] == unique(c_s_df[,strata])[i], cluster]
          clust_samp <- c(boot_clusts, sample(temp_clusts, size = length(temp_clusts) - 1, replace = TRUE))
          # length(temp_clusts) - 1 because setting m_h = n_h - 1 means we don't need
          # a weight correction
        }
        
        boot_clusts <- clust_samp
        
        # temporary boot_df with a single observation for each sampled cluster
        df_boot <- df[df[,cluster] %in% boot_clusts,]
        
        # generate bootstrap dataframe using these clusters
        for (k in 1:length(unique(boot_clusts))) {
          which_clust <- unique(boot_clusts)[k]
          
          # how many times is this cluster repeated
          num_repeated <- length(which(boot_clusts == which_clust))
          
          # if num_repeated > 1, then add in those repeats to df_boot
          if (num_repeated > 1) {
            temp_df <- df[df[,cluster] == which_clust,]
            df_boot <- rbind(df_boot,temp_df[rep(1:nrow(temp_df),num_repeated - 1),])
          }
        }
        
        # add correction factor to design weights
        # NO CORRECTION NEEDED if sample m_h = n_h - 1 from each strata h
        # df_boot <- suppressMessages(left_join(df_boot, h_df)) 
        # df_boot$correction <- df_boot$n_h / (df_boot$n_h - 1)
        # df_boot[,weights] <- df_boot[,weights] * df_boot$correction
        
        # get turnbull estimator on bootstrapped dataset
        small_samp_est <- pssst:::rcpp_turnbull(niter = niter, 
                                                t0 = df_boot[,t_0i],
                                                t1 = df_boot[,t_1i],
                                                lefttrunc = df_boot[,lefttrunc], 
                                                righttrunc = df_boot[,righttrunc],
                                                weights = df_boot[,weights],
                                                set_lower = res_full[[1]],
                                                set_upper = res_full[[2]])
        
        # fill in boot_mat
        boot_mat[,i] <- small_samp_est[[3]][,niter]
      }
      
      # calculate bootstrap variance
      boot_var <- apply(boot_mat,1,var)
    }
    
    # create final data frame
    if (bootstrap_var) {
      res_df <- data.frame(t0 = res[[1]],
                           t1 = res[[2]],
                           est = 1 - cumsum(res[[3]][,niter]),
                           var_est = boot_var)
      ret_lst <- list(result = res_df,
                      iterations = apply(res[[3]],2,function(x) {1 - cumsum(x)}),
                      bootstrap_samps = apply(boot_mat,2,function(x) {1 - cumsum(x)}),
                      niter = niter,
                      niter_bootstrap = niter_bootstrap)
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
    period_lst <- vector("list", length(unique(df$period)))
    
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
      
      # calculate bootstrap variance
      if (bootstrap_var) {
        message(paste0("Calculating bootstrap variance for period ",sort(unique(df$period))[l],"..."))
        
        # get data frame containing unique clusters and stratas
        c_s_df <- df_p[,c(cluster, strata)] %>% dplyr::distinct()
        
        # number of strata
        h_df <- df_p %>% 
          dplyr::group_by(strata) %>% 
          dplyr::summarise(n_h = length(unique(cluster))) %>% 
          as.data.frame()
        h <- nrow(h_df)
        
        # vector of unique clusters
        unique_clusts <- unique(df_p[,cluster])
        
        # set up matrix for bootstrap results
        boot_mat <- matrix(NA, nrow = length(res_full[[1]]), ncol = niter_bootstrap)
        
        # for i in 1:niter_bootstrap...
        # - sample clusters with replacement
        # - multiply design weights by correction factor n_h / (n_h - 1)
        # - calculate turnbull estimator
        
        for (i in 1:niter_bootstrap) {
          # sample clusters with replacement
          boot_clusts <- sample(unique_clusts, size = length(unique_clusts), replace = TRUE)
          
          # temporary boot_df with a single observation for each sampled cluster
          df_boot <- df_p[df_p[,cluster] %in% boot_clusts,]
          
          # generate bootstrap dataframe using these clusters
          for (k in 1:length(unique(boot_clusts))) {
            which_clust <- unique(boot_clusts)[k]
            
            # how many times is this cluster repeated
            num_repeated <- length(which(boot_clusts == which_clust))
            
            # if num_repeated > 1, then add in those repeats to df_boot
            if (num_repeated > 1) {
              temp_df <- df_p[df_p[,cluster] == which_clust,]
              df_boot <- rbind(df_boot,temp_df[rep(1:nrow(temp_df),num_repeated - 1),])
            }
          }
          
          # add correction factor to design weights
          df_boot <- suppressMessages(left_join(df_boot, h_df)) 
          df_boot$correction <- df_boot$n_h / (df_boot$n_h - 1)
          df_boot[,weights] <- df_boot[,weights] * df_boot$correction
          
          # get turnbull estimator on bootstrapped dataset
          small_samp_est <- pssst:::rcpp_turnbull(niter = niter, 
                                          t0 = df_boot[,t_0i],
                                          t1 = df_boot[,t_1i],
                                          lefttrunc = df_boot[,lefttrunc], 
                                          righttrunc = df_boot[,righttrunc],
                                          weights = df_boot[,weights],
                                          set_lower = res_full[[1]],
                                          set_upper = res_full[[2]])
          
          # fill in boot_mat
          boot_mat[,i] <- small_samp_est[[3]][,niter]
        }
        
        # calculate bootstrap variance
        boot_var <- apply(boot_mat,1,var)
        
      }
      
      # create final data frame
      if (bootstrap_var) {
        res_df <- data.frame(t0 = res[[1]],
                             t1 = res[[2]],
                             est = 1 - cumsum(res[[3]][,niter]),
                             var_est = boot_var)
        ret_lst <- list(result = res_df,
                        iterations = apply(res[[3]],2,function(x) {1 - cumsum(x)}),
                        bootstrap_samps = apply(boot_mat,2,function(x) {1 - cumsum(x)}),
                        niter = niter,
                        niter_bootstrap = niter_bootstrap)
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
    if (bootstrap_var) {
      
      # combine results from all periods into one long dataframe
      temp_lst <- vector("list", length = length(period_lst))
      temp_lst_iter <- vector("list", length = length(period_lst))
      temp_lst_boot <- vector("list", length = length(period_lst))
      for (i in 1:length(period_lst)) {
        period_lst[[i]][[1]][,period] <- sort(unique(df$period))[i]
        temp_lst[[i]] <- period_lst[[i]][[1]][,c(period, "t0", "t1", "est", "var_est")]
        temp_lst_iter[[i]] <- period_lst[[i]][[2]]
        temp_lst_boot[[i]] <- period_lst[[i]]$bootstrap_samps
      }
      res_df <- do.call(rbind, temp_lst)
      
      ret_lst <- list(result = res_df,
                      iterations = temp_lst_iter,
                      bootstrap_samps = temp_lst_boot,
                      niter = niter,
                      niter_bootstrap = niter_bootstrap)
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

