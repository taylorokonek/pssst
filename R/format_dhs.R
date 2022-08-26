#' Benchmark posterior draws via a Metropolis-Hastings algorithm
#' 
#' Benchmarks draws from a posterior distribution using a Metropolis-Hastings approach. Takes in 
#' a matrix of posterior draws, a matrix of fixed effect draws, and national-level benchmark data.
#' 
#' @param df a dataframe containing a births recode file from DHS, read into R beforehand using
#' \code{readstata13::read.dta13}
#' @param survey_year the (numeric) year of the DHS survey. If survey spand multiple years,
#' enter the latest one
#' @param period_boundaries a numeric vector containing the years surrounding each time period of interest.
#' Example: if periods are 2000-2005 and 2005-2010, \code{period_boundaries} should be \code{c(2000,2005,2010)}.
#' @param strata a string vector containing which column in \code{df} contains the strata 
#' information. Can be multiple columns. Defaults to "v023".
#' @return A list containing: 
#' \itemize{
#' \item fitted_list: a list of matrices of benchmarked posterior samples of fitted values in 
#' order arrange(time). Each matrix will have rows arranged in order arrange(region).
#' \item natl_list: a list of vectors containing aggregated national-level samples that were 
#' accepted, in order arrange(time)
#' \item prop_accepted: the proportion of samples accepted during sampling. 
#' } 
#' 
#' @author Taylor Okonek
#' @export format_dhs

format_dhs <- function(df, 
                       survey_year, 
                       period_boundaries, 
                       strata = c("v023", "v024", "v025")[1]) {
  
  # get year_cut from period_boundaries
  year_cut <- period_boundaries[-length(period_boundaries)]
  
  # call get_births
  births <- get_births(dat = df, surveyyear = survey_year, year.cut = year_cut)
  
  # remove individuals born before min(year_cut)
  num_before <- length(which(births$year_born < min(year_cut)))
  message(paste0("Removing ", num_before, " children born before ", min(year_cut)))
  births <- births %>% filter(year_born >= min(year_cut))
  
  # expand dataframe
  births <- df_expand(surv_df = births, 
                      period_boundaries = period_boundaries, 
                      survey_date = survey_year)
  
  return(births)
}