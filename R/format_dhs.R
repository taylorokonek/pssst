#' Formats DHS survey data in appropriate form for parametric survival model for synthetic children.
#' 
#' Formats DHS survey data in appropriate way for input to \code{surv_synthetic}. Takes in 
#' a dataframe containing a births recode file from DHS.
#' 
#' @param df a dataframe containing a births recode file from DHS, read into R beforehand using
#' \code{readstata13::read.dta13}
#' @param survey_year the (numeric) year of the DHS survey. If survey spand multiple years,
#' enter the latest one
#' @param period_boundaries a numeric vector containing the years surrounding each time period of interest.
#' Example: if periods are 2000-2005 and 2005-2010, \code{period_boundaries} should be \code{c(2000,2005,2010)}.
#' @param strata a string vector containing which column in \code{df} contains the strata 
#' information. Can be multiple columns. Defaults to "v023".
#' @return A dataframe containing the births recode in a format that can be input to \code{surv_synthetic}.
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
  
  # remove people who are right-censored at age 0 (absolute nonsense)
  num_nonsense <- births %>% filter(I_i == 0 & t_i == 0) %>% nrow()
  message(paste0("Removing ", num_nonsense, " children right-censored at age 0 "))
  births <- births %>% filter(!(I_i == 0 & t_i == 0))
  
  # divide weights by constant because DHS
  births <- births %>% mutate(weights = weights/1e6)
  
  return(births)
}