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
#' @param right_censor_time age in months at which every individual who was not yet died
#' should be right-censored. If there is no such age, set equal to NA. If estimating U5MR, 
#' consider setting to 60.
#' @param intervals an optional character vector specifying any specific intervals in which
#' observations should be censored. For example, to interval censor people between ages 6 and 18
#' months, \code{intervals} = c("0-1","1-2,"6-18"). Intervals cannot overlap (i.e. can't have 6-18, 12-20). Default
#' intervals are the ones observed in DHS: exact daily deaths before 1 month, monthly through age 24 months, yearly after.
#' Days will be reported as 1/30th of a month. If an individual "Died on day of birth", they are
#' interval censored from [0, 1/30]
#' @param strata a string vector containing which column in \code{df} contains the strata 
#' information. Can be multiple columns. Defaults to "v023".
#' @return A dataframe containing the births recode in a format that can be input to \code{surv_synthetic}. Each
#' individual will have multiple rows in the dataframe (one for each period).
#' \itemize{
#'  \item individual: individual id
#'  \item household: household id
#'  \item cluster: cluster id
#'  \item strata
#'  \item weights: survey weights
#'  \item p: numeric period identified
#'  \item y_p: the date the period begins
#'  \item l_p: the length of the period, in months
#'  \item b_i: the birth date of the individual
#'  \item t_i: age at right-censoring, if right-censored
#'  \item t_0i: lower bound of interval, if interval-censored. date of death if observed exactly
#'  \item t_1i: upper bound of interval, if interval-censored. date of death if observed exactly
#'  \item I_i: indicator for death. I_i == 1 if interval-censored or died exactly, I_i == 0 if right censored
#'  \item survey_year: year of survey
#'  \item a_pi: age of child in months at beginning of time period
#'  \item exact_age: exact age of death recorded in DHS (b6 in births recode)
#'  \item A_i: indicator for if a child is interval-censored across the boundary of a time period
#' }
#' 
#' @author Taylor Okonek
#' @export format_dhs
format_dhs <- function(df, 
                       survey_year, 
                       period_boundaries, 
                       right_censor_time = NA,
                       intervals = NA,
                       strata = c("v023", "v024", "v025")[1]) {
  
  max_year <-period_boundaries[length(period_boundaries)]
  
  # get year_cut from period_boundaries
  year_cut <- period_boundaries[-length(period_boundaries)]
  
  # call get_births
  births <- get_births(dat = df, 
                       surveyyear = survey_year, 
                       year.cut = year_cut, 
                       strata = strata,
                       intervals = intervals)
  
  # remove individuals who died before min(year_cut)
  num_before <- length(which(births$year_died < min(year_cut)))
  
  if(num_before > 0) {
    message(paste0("Removing ", num_before, " children who died before ", min(year_cut)))
  }
  births <- births %>% filter(year_died >= min(year_cut))
  
  # remove individuals born after max_year
  num_after <- length(which(births$year_born >= max_year))
  if (num_after > 0) {
    message(paste0("Removing ", num_after, " children born after ", max_year))
  }
  births <- births %>% filter(year_born < max_year)
  
  # if a right_censor_time is specified, do that
  if (!is.na(right_censor_time)) {
    # identify who dies after right_censor_time
    row_ids <- which(births$died & (births$age_at_censoring > right_censor_time))
    
    if (length(row_ids) > 0) {
      birth_dates <- paste(births[row_ids,]$year_born, 
                           births[row_ids,]$month_born, sep = "-") %>% ym()
      right_cens_dates <- birth_dates %m+% months(right_censor_time)
      
      # identify who would be right censored before the first period
      rc_before_ids <- which(right_cens_dates < (paste(min(year_cut),"1",sep = "-") %>% ym() %m+% months(1)))
      rc_before_ids_final <- row_ids[rc_before_ids]
      
      # instead make them right censored
      births[row_ids,]$right_censored <- TRUE
      births[row_ids,]$interval_censored <- FALSE
      births[row_ids,]$died <- FALSE
      births[row_ids,]$age_at_censoring <- right_censor_time
      births[row_ids,]$t0 <- 0
      births[row_ids,]$t1 <- Inf
    }
    
    # identify who is right censored later than right_censor_time
    row_ids <- which(!births$died & (births$age_at_censoring > right_censor_time))
    
    if (length(row_ids) > 0) {
      births[row_ids,]$age_at_censoring <- right_censor_time
    }
    
    # remove individuals who are right-censored before first time period based on cutoff
    if (length(rc_before_ids_final) > 0) {
      message(paste0("Removing ", length(rc_before_ids_final), " children who are right-censored before ", min(year_cut)))
      births <- births[-rc_before_ids_final,]
    }
    
  }
  
  
  # death with exact days of death
  # get rows where child's death is recorded in days
  suppressWarnings(temp_df <- births %>% 
                     dplyr::mutate(exact_age = as.numeric(as.character(exact_age))))
  exact_rows <- which(temp_df$exact_age < 200 & temp_df$exact_age >= 100)
  exact_rows <- c(exact_rows,which(births$exact_age == "Days: 1"))
  
  # get birth in days
  daily_births <- births[exact_rows,]$exact_age 
  suppressWarnings(daily_births <- ifelse(daily_births == "Days: 1", 1, as.numeric(as.character(daily_births)) - 100))
  daily_births <- daily_births / 30 # transform to months
  
  # if any exact births are within intervals specified, don't alter them
  which_remove <- c()
  lbs <- c()
  for (i in 1:length(intervals)) {
    lb <- intervals[i] %>% str_split("-") %>% unlist %>% nth(1) %>% as.numeric
    ub <- intervals[i] %>% str_split("-") %>% unlist %>% nth(2) %>% as.numeric
    which_remove <- c(which_remove, which((daily_births <= ub) & (daily_births >= lb)))
    lbs <- c(lbs, lb)
  }
  which_remove <- unique(which_remove)
  
  if (length(which_remove) > 0) {
    exact_rows <- exact_rows[-which_remove]
  }
  
  # replace t0 and t1 in dataframe with daily_births
  births[exact_rows,]$t0 <- daily_births
  births[exact_rows,]$t1 <- daily_births
  
  # censor children who "Died on day of birth"
  if (!(0 %in% lbs)) {
    num_died_at_birth <- births %>% filter(exact_age == "Died on day of birth") %>% nrow()
    message("Interval censoring ",num_died_at_birth, " children who died on day of birth from [0,1] day")
    which_died <- which(births$exact_age == "Died on day of birth")
    births[which_died,]$t1 <- 1/30
    births[which_died,]$age_at_censoring <- 1/30
  }
  
  # expand dataframe
  births <- df_expand(surv_df = births, 
                      period_boundaries = period_boundaries, 
                      survey_date = survey_year)
  
  # if a death is observed exactly, t_0i == t_1i, and I_i == 1
  
  # remove people who are right-censored at age 0 (absolute nonsense)
  num_nonsense <- births %>% filter(I_i == 0 & t_i == 0) %>% nrow()
  if (num_nonsense > 0) {
    message(paste0("Removing ", num_nonsense, " children right-censored at age 0 "))
  }
  births <- births %>% filter(!(I_i == 0 & t_i == 0))
  
  # divide weights by constant because DHS
  births <- births %>% mutate(weights = weights/1e6)
  
  return(births)
}