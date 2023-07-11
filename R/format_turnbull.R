#' Formats output from \code{format_dhs} in appropriate form for \code{turnbull} 
#' 
#' Formats output from \code{format_dhs} in appropriate way for input to \code{turnbull}. Takes in 
#' the exact output of \code{format_dhs}
#' 
#' @param df output of \code{format_dhs}
#' @return A dataframe containing the births data in a format that can be input to \code{turnbull}.
#' 
#' @author Taylor Okonek
#' @export format_turnbull
format_turnbull <- function(df) {
  
  # All of this should be moved to c++ because it is VERY slow
  
  # get data ready for turnbull function
  I0_df <- df %>% dplyr::filter(I_i == 0)
  I1_df <- df %>% dplyr::filter(I_i == 1)
  
  I0_df$lefttrunc <- NA
  I0_df$righttrunc <- NA
  
  I1_df$lefttrunc <- NA
  I1_df$righttrunc <- NA
  
  # set up dataframe
  ret_small_0 <- I0_df[1,]
  
  for (i in 1:length(unique(I0_df$individual))) {
    row_ids <- which(I0_df$individual == unique(I0_df$individual)[i])
    test <- I0_df[row_ids,]
    test <- test %>%
      dplyr::mutate(lefttrunc = ifelse(a_pi < 0, 0, a_pi),
             righttrunc = Inf,
             t_0i = lefttrunc)
    test$t_1i <- c(test$lefttrunc[-1],test$t_i[1])
    
    # get rows where they are alive
    test_small <- test[(test$a_pi >= -test$l_p) & (test$a_pi < test$t_i[1]),]
    
    # for turnbull() function, t_0i needs to be censoring time, t_1i needs to be Inf 
    test_small$t_1i <- ifelse(test_small$t_1i > test_small$t_i[1], test_small$t_i[1], test_small$t_1i)
    test_small$t_0i <- test_small$t_1i
    test_small$t_1i <- Inf
    
    # rbind
    ret_small_0 <- rbind(ret_small_0, test_small)
  }
  
  # set up dataframe
  ret_small <- I1_df[1,]
  
  for (i in 1:length(unique(I1_df$individual))) {
    row_ids <- which(I1_df$individual == unique(I1_df$individual)[i])
    test <- I1_df[row_ids,]
    
    A_i_orig <- test$A_i[1]
    
    # if not interval censored across time period boundary
    if (A_i_orig == 0) {
      # which time period are they interval censored in
      which_cens <- tail(which(test$t_0i >= test$a_pi), 1)
      
      # which time periods are they alive in
      which_alive <- which(((test$a_pi + test$l_p) > 0) & !((test$t_1i + test$l_p) <= (test$a_pi + test$l_p)))
      
      # if they're only alive in one period...
      if (length(which_alive == 1)) {
        test_small <- test[which_cens,]
        test_small$lefttrunc <- 0
        test_small$righttrunc <- Inf
        
        # rbind
        ret_small <- rbind(ret_small, test_small)
      } else {
        test_small <- test[which_cens,]
        test_small$lefttrunc <- pmax(0, test_small$a_pi)
        test_small$righttrunc <- Inf
        
        test_small_right <- test[which_alive[!(which_alive %in% which_cens)],]
        test_small_right$t_0i <- test_small_right$a_pi + test_small_right$l_p
        test_small_right$t_1i <- Inf
        test_small_right$I_i <- 0
        test_small_right$A_i <- 0
        
        test_small_right$lefttrunc <- pmax(test_small_right$a_pi,0)
        test_small_right$righttrunc <- Inf
        
        test_small <- rbind(test_small_right, test_small)
        
        # rbind
        ret_small <- rbind(ret_small, test_small)
      }
    
    # if interval censored across time period boundary
    } else {
      # which time periods are they interval censored in
      which_cens <- which((test$t_0i <= (test$a_pi + test$l_p)) & (test$t_1i > test$a_pi))
      
      if (length(which_cens) != 2) {
        stop("function not built to handle people interval censored across three time periods yet")
      }
      
      # which time periods are they alive in
      which_alive <- which(((test$a_pi + test$l_p) > 0) & !((test$t_1i + test$l_p) <= (test$a_pi + test$l_p)))
      
      # get the interval censored periods
      test_small <- test[which_cens,]
      test_small$t_0i <- c(test_small$t_0i[1], test_small$a_pi[2])
      test_small$t_1i <- c(test_small$a_pi[2], test_small$t_1i[1])
      test_small$lefttrunc <- pmax(0, test_small$a_pi)
      test_small$righttrunc <- Inf
      
      # re-weight A_i = 1 intervals by length
      interval_lengths <- test_small$t_1i - test_small$t_0i
      test_small$weights <- test_small$weights * interval_lengths/sum(interval_lengths)
      
      # add in their right-censored times for earlier periods (if they have any)
      if (length(which_alive) > 2) {
        test_small_right <- test[which_alive[!(which_alive %in% which_cens)],]
        test_small_right$t_0i <- test_small_right$a_pi + test_small_right$l_p
        test_small_right$t_1i <- Inf
        test_small_right$I_i <- 0
        test_small_right$A_i <- 0
        
        test_small_right$lefttrunc <- pmax(test_small_right$a_pi,0)
        test_small_right$righttrunc <- Inf
        
        test_small <- rbind(test_small_right, test_small)
      }
      
      # rbind
      ret_small <- rbind(ret_small, test_small)
    }
    

  }
  
  # remove first row of ret_small
  ret_small_0 <- ret_small_0[-1,]
  ret_small <- ret_small[-1,]
  
  turnbull_df <- rbind(ret_small_0, ret_small)
  
  return(turnbull_df)
}