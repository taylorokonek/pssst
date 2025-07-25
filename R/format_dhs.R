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
#' months, \code{intervals} = c("0-1","1-2","6-18"). Intervals cannot overlap (i.e. can't have 6-18, 12-20). Default
#' intervals are the ones observed in DHS: exact daily deaths before 1 month, monthly through age 24 months, yearly after.
#' Days will be reported as 1/30th of a month. If an individual "Died on day of birth", they are
#' interval censored from [0, 1/30]
#' @param cmc_adjust number of months to add to the recorded month in the dataset. As an example, 
#' the Ethiopian calendar is 92 months behind the Gregorian calendar in general, so if a DHS survey
#' from Ethiopia is used, \code{cmc_adjust} should be set to 92. Default value is 0.
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
                       cmc_adjust = 0,
                       strata = c("v023", "v024", "v025")[1]) {
  
  max_year <- period_boundaries[length(period_boundaries)]
  
  # get year_cut from period_boundaries
  year_cut <- period_boundaries[-length(period_boundaries)]
  
  # if any of the necessary variables are haven_labelled, change them
  if (class(df$b5)[1] == "haven_labelled") {
    # convert some variables to factors
    alive <- attr(df$b5, which = "labels")
    names(alive) <- tolower(names(alive))
    df$b5 <- ifelse(unclass(df$b5) == alive["yes"][[1]], "yes", "no")
    df$b5 <- factor(df$b5, levels = c("yes", "no"))
  }
  
  if (class(df$b6)[1] == "haven_labelled") {
    # convert some variables to factors
    exact_labels <- attr(df$b6, which = "labels")
    names(exact_labels) <- tolower(names(exact_labels))
    
    # the only two we care about are died on day of birt
    temp <- unclass(df$b6)
    for (j in 1:length(exact_labels)) {
      temp <- ifelse(temp == exact_labels[j], names(exact_labels)[j], temp)
    }
    df$b6 <- temp
  }
  
  if ("v025" %in% strata) {
    if (class(df$v025)[1] == "haven_labelled") {
      strat <- attr(df$v025,which='labels')
      names(strat) <- tolower(names(strat))
      df$v025 <- ifelse(unclass(df$v025) == strat["urban"][[1]],'urban','rural')
      df$v025 <- factor(df$v025, levels = c('urban','rural'))
    }
  }
  
  if ("v023" %in% strata) {
    if (class(df$v023)[1] == "haven_labelled") {
      # check if there are any strata with no data
      strata_names <- attr(df$v023, which = "labels") %>% names()
      strata_names <- strata_names[strata_names != "missing"]
      
      # all possible strata values
      strata_vals <- attr(df$v023, which = "labels") %>% unname()
      
      # strata with observations in data frame
      df$v023 <- df$v023 %>% unclass()
      obs_vals <- df$v023 %>% unique %>% sort()
      which_missing <- strata_vals[!(strata_vals %in% obs_vals)]
      
      if (length(which_missing) > 0) {
        df$v023 <- factor(df$v023, levels = df$v023 %>% table() %>% names(),
                          labels = strata_names[-which_missing])
      } else {
        df$v023 <- factor(df$v023, levels = df$v023 %>% table() %>% names(),
                          labels = strata_names)
      }
    }
  }
  
  if ("v022" %in% strata) {
    if (class(df$v022)[1] == "haven_labelled") {
      df$v022 <- df$v022 %>% unclass()
      df$v022 <- factor(df$v022, levels = df$v022 %>% table() %>% names(),
                        labels = attr(df$v022, which = "labels") %>% names())
    }
  }
  
  if ("v024" %in% strata) {
    if (class(df$v024)[1] == "haven_labelled") {
      df$v024 <- df$v024 %>% unclass()
      df$v024 <- factor(df$v024, levels = df$v024 %>% table() %>% names(),
                        labels = attr(df$v024, which = "labels") %>% names())
    }
  }
  
  # call get_births
  births <- get_births(dat = df, 
                       surveyyear = survey_year, 
                       year.cut = year_cut, 
                       strata = strata,
                       cmc.adjust = cmc_adjust,
                       intervals = intervals)
  
  # convert things to vectors (if needed, as would be the case if data comes from rdhs)
  for (i in 1:ncol(births)) {
    births[,i] <- as.vector(births[,i])
  }
  
  # remove individuals who died before min(year_cut)
  num_before <- length(which(births$year_died < min(year_cut)))
  if(num_before > 0) {
    message(paste0("Removing ", num_before, " children who died before ", min(year_cut)))
  }
  births <- births[which(is.na(births$year_died) | (births$year_died >= min(year_cut))),]
  
  # remove individuals who are right-censored before min(year_cut)
  birth_dates <- suppressWarnings((paste(births$year_born, births$month_born, sep = "-") %>% ym()))
  right_censored_dates <- birth_dates %m+% months(births$age_at_censoring)
  censored_before <- which(right_censored_dates <= ym(paste0(min(year_cut), "-01")))
  if (length(censored_before) > 0) {
    message(paste0("Removing ", length(censored_before), " children who were right censored before ", min(year_cut)))
    births <- births[-censored_before,]
  }
  
  # if right_censor_time specified and they reached this age before min(year), remove them
  if (!is.na(right_censor_time)) {
    right_censored_dates_new <- birth_dates %m+% months(right_censor_time)
    censored_before_new <- which(right_censored_dates_new <= ym(paste0(min(year_cut), "-01")))
    if (length(censored_before_new) > 0) {
      message(paste0("Removing ", length(censored_before_new), " children who would be right censored at age ",right_censor_time," before ", min(year_cut)))
      births <- births[-censored_before_new,]
    }
  }
  
  # remove individuals born after max_year
  num_after <- length(which(births$year_born >= max_year))
  if (num_after > 0) {
    message(paste0("Removing ", num_after, " children born after ", max_year))
  }
  births <- births %>% filter(year_born < max_year)
  
  # if they died after max_year, right censor them at their age at max year
  died_dates <- suppressWarnings((paste(births$year_died, births$month_died, sep = "-") %>% ym()))
  max_date <- ym(paste0(max_year,"-01"))
  died_after_max <- which((!is.na(died_dates)) & (died_dates >= max_date))
  
  new_rightcensoringage <- floor(-((ym(paste(births$year_born, births$month_born, sep = "-")) - 
                                      ym(paste0(max_year,"-01")))[died_after_max]) / 30)
  
  row_ids <- died_after_max
  # instead make them right censored
  if (length(row_ids) > 0) {
    births[row_ids,]$right_censored <- TRUE
    births[row_ids,]$interval_censored <- FALSE
    births[row_ids,]$died <- FALSE
    births[row_ids,]$age_at_censoring <- new_rightcensoringage
    births[row_ids,]$t0 <- 0
    births[row_ids,]$t1 <- Inf
  }
  
  # if a right_censor_time is specified, do that
  if (!is.na(right_censor_time)) {
    # identify who dies after right_censor_time
    row_ids1 <- which(births$died & (births$age_at_censoring > right_censor_time))
    
    rc_before_ids_final <- c()
    rc_before_ids_final2 <- c()
    
    if (length(row_ids1) > 0) {
      birth_dates <- paste(births[row_ids1,]$year_born, 
                           births[row_ids1,]$month_born, sep = "-") %>% ym()
      right_cens_dates <- birth_dates %m+% months(right_censor_time)
      
      # identify who would be right censored before the first period
      rc_before_ids <- which(right_cens_dates < (paste(min(year_cut),"1",sep = "-") %>% ym() %m+% months(1)))
      rc_before_ids_final <- row_ids1[rc_before_ids]
      
      # instead make them right censored
      births[row_ids1,]$right_censored <- TRUE
      births[row_ids1,]$interval_censored <- FALSE
      births[row_ids1,]$died <- FALSE
      births[row_ids1,]$age_at_censoring <- right_censor_time
      births[row_ids1,]$t0 <- 0
      births[row_ids1,]$t1 <- Inf
    }
    
    # identify who is right censored later than right_censor_time
    row_ids2 <- which(!births$died & (births$age_at_censoring > right_censor_time))
    
    if (length(row_ids2) > 0) {
      births[row_ids2,]$age_at_censoring <- right_censor_time
      
      birth_dates2 <- paste(births[row_ids2,]$year_born, 
                            births[row_ids2,]$month_born, sep = "-") %>% ym()
      right_cens_dates2 <- birth_dates2 %m+% months(right_censor_time)
      
      # identify who would be right censored before the first period
      rc_before_ids2 <- which(right_cens_dates2 < (paste(min(year_cut),"1",sep = "-") %>% ym() %m+% months(1)))
      rc_before_ids_final2 <- row_ids2[rc_before_ids2]
    }
    
    rc_before_ids_final <- unique(c(rc_before_ids_final, rc_before_ids_final2))
    rc_before_ids_final <- rc_before_ids_final[!is.na(rc_before_ids_final)]
    
    # remove individuals who are right-censored before first time period based on cutoff
    if (length(rc_before_ids_final) > 0) {
      message(paste0("Removing ", length(rc_before_ids_final), " children who are right-censored before ", min(year_cut)))
      births <- births[-rc_before_ids_final,]
    }
    
  }
  
  # deal with special case (exact age 199 set to 0 months)
  births$exact_age <- ifelse(births$exact_age == 199, 200, births$exact_age)
  
  # death with exact days of death
  # get rows where child's death is recorded in days
  suppressWarnings(temp_df <- births %>% 
                     dplyr::mutate(exact_age = as.numeric(as.character(exact_age))))
  exact_rows <- which(temp_df$exact_age < 200 & temp_df$exact_age >= 100)
  exact_rows <- c(exact_rows, which(births$exact_age %in% c("Days: 1", "days: 1")))
  
  # get birth in days
  daily_births <- births[exact_rows,]$exact_age 
  suppressWarnings(daily_births <- ifelse(daily_births %in% c("Days: 1", "days: 1"), 1, as.numeric(as.character(daily_births)) - 100))
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
    daily_births <- daily_births[-which_remove]
  }
  
  # replace t0 and t1 in dataframe with daily_births
  if (length(exact_rows) > 0) {
    births[exact_rows,]$t0 <- daily_births
    births[exact_rows,]$t1 <- daily_births
  }
  
  # censor children who "Died on day of birth"
  if (!(0 %in% lbs)) {
    num_died_at_birth <- births %>% filter(exact_age == "Died on day of birth" |
                                             exact_age == "died on day of birth" |
                                             exact_age == 100) %>% nrow()
    if (num_died_at_birth > 0) {
      message("Interval censoring ",num_died_at_birth, " children who died on day of birth from [0,1] day")
      which_died <- which(births$exact_age == "Died on day of birth" |
                            births$exact_age == "died on day of birth" |
                            births$exact_age == 100)
      births[which_died,]$t1 <- 1/30
      births[which_died,]$age_at_censoring <- 1/30
    }
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
  
  # deal with people interval censored across last boundary
  uniq_indivs <- unique(births$individual)
  change_As <- c()
  
  for (i in 1:length(uniq_indivs)) {
    temp <- births[which(births$individual == uniq_indivs[i]),]
    
    # if right-censored, make sure t_i <= max(a_pi + l_p)
    if (temp$I_i[1] == 0) {
      new_t_i <- min(temp$t_i[1],max(temp$a_pi + temp$l_p))
      births[which(births$individual == uniq_indivs[i]),]$t_i <- new_t_i
      
      # if interval-censored or exactly observed, make sure 
    } else {
      new_t_1i <- min(temp$t_1i[1],max(temp$a_pi + temp$l_p))
      new_t_0i <- max(temp$t_0i[1], min(temp$a_pi))
      
      if ((new_t_1i != temp$t_1i[1]) | (new_t_0i != temp$t_0i[1])) {
        change_As <- c(change_As, which(births$individual == uniq_indivs[i]))
      }
      
      births[which(births$individual == uniq_indivs[i]),]$t_1i <- new_t_1i
      births[which(births$individual == uniq_indivs[i]),]$t_0i <- new_t_0i
    }
  }
  
  if(length(change_As) > 0) {
    births[change_As,]$A_i <- 0
  }
  
  # make sure I_i is exactly 1 or 0
  births$I_i <- round(births$I_i)
  
  # make sure household, individual, cluster, and strata are all integers
  births$individual <- round(births$individual)
  births$household <- round(births$household)
  births$cluster <- round(births$cluster)
  #births$strata <- round(as.numeric(births$strata))
  
  return(births)
}