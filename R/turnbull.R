#' Function to obtain Turbull estimates for arbitrarily censored and truncated observations
#' 
#' @param df a dataframe containing the output from \code{format_dhs}, or optionally, dataframe
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
#' @param niter number of iterations to run the algorithm
#' @return a list containing the following:
#' \begin{itemize}
#' \item a dataframe of final estimates for specific times
#' \item a matrix of estimates at each iteration, with ncol = niter
#' \item number of iterations
#' \end{itemize}
#' 
#' @author Taylor Okonek
#' @export turnbull
turnbull <- function(df, 
                     t_0i = "t_0i",
                     t_1i = "t_1i", 
                     lefttrunc = "lefttrunc", 
                     righttrunc = "righttrunc", 
                     weights = "weights", 
                     niter) {
  # error checking - TO DO
  
  # call rcpp_turnbull
  res <- rcpp_turnbull(niter = niter, 
                       t0 = df[,t_0i], 
                       t1 = df[,t_1i], 
                       lefttrunc = df[,lefttrunc], 
                       righttrunc = df[,righttrunc], 
                       weights = df[,weights]) 
  
  # create final data frame
  res_df <- data.frame(t0 = res[[1]],
             t1 = res[[2]],
             est = 1 - cumsum(res[[3]][,niter]))
  
  return(list(result = res_df,
              iterations = apply(res[[3]],2,function(x) {1 - cumsum(x)}),
              niter = niter))
}

