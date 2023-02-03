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
  
  for (i in 1:length(unique(I0_df$individual))) {
    row_ids <- which(I0_df$individual == I0_df$individual[i])
    test <- I0_df[row_ids,]
    test <- test %>%
      mutate(lefttrunc = ifelse(a_pi < 0, 0, a_pi),
             righttrunc = Inf,
             t_0i = lefttrunc)
    test$t_1i <- c(test$lefttrunc[-1],test$t_i[1])
    I0_df[row_ids,] <- test
  }
  
  I0_df_small <- I0_df %>% dplyr::filter(a_pi >= -60) 
  
  for (i in 1:length(unique(I1_df$individual))) {
    row_ids <- which(I1_df$individual == unique(I1_df$individual)[i])
    test <- I1_df[row_ids,]
    test <- test %>%
      mutate(lefttrunc = ifelse(a_pi < 0, 0, a_pi),
             righttrunc = Inf,
             t_0i = pmin(t_0i,lefttrunc))
    test$t_1i <- pmin(c(test$lefttrunc[-1],test$t_1i[1]), test$t_1i[1])
    I1_df[row_ids,] <- test
  }
  
  I1_df_small <- I1_df %>% 
    filter(lefttrunc < t_1i)
  
  turnbull_df <- rbind(I0_df_small, I1_df_small)
  
  return(turnbull_df)
}