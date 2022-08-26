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
#' @return TBD
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
                           t_1i = "t_1i") {
  
  # make new df with appropriate columns
  temp <- df[,c(individual, household, cluster, strata, weights, p, a_pi, l_p,
                I_i, A_i, t_i, t_0i, t_1i)]
  df <- temp
  
  # pivot wider
  df <- df %>%
    pivot_wider(id_cols = c(individual, household, cluster, strata, weights, I_i, A_i, t_i, t_0i, t_1i),
                names_from = p,
                values_from = c(a_pi, l_p))

  return(df)
}