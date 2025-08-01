#' Function to get cleaned births data frame from a DHS births dataframe
#' 
#' @param dat a births data frame
#' @param surveyyear numeric, year of survey
#' @param strata string vector, which columns in dat correspond to strata
#' @param year.cut a numeric vector containing the first year of each time period of interest
#' @return organized births dataframe
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
get_births <- function(dat, 
                       surveyyear, 
                       strata, 
                       year.cut,
                       cmc.adjust = 0,
                       intervals) {
  
  # additional definitions
  variables <- c("caseid", 
                "v001", "v002", "v004", "v005", "v021", "v022", "v023", 
                "v024", "v025", "v139", "bidx")
  dob <- "b3"
  alive <- "b5"
  age <- "b7"
  exact_age <- "b6"
  age.truncate <- 24
  # cmc.adjust <- 0
  date.interview <- "v008"
  month.cut <- c(seq(1,24), seq(36,12*100, by = 12))
  
  # DHS I is missing v021 (psu), v022 (strata), v023 (other strata variable), v024 (de facto region of residence),
  # v025 (de facto type of place of residence, i.e. urban/rural), v139 (de jure region of usual residence)
  # -- these columns aren't critical, so just add as NA
  missing_cols <- setdiff(c("v021", "v022", "v023", "v024", "v025", "v139"), names(dat))
  dat[missing_cols] <- NA
  
  # if strata is NULL or NA, add empty column
  if (is.null(strata) || all(is.na(strata))) {
    dat$strata <- NA
    strata <- "strata"
  }
  
  # extra interval censoring, if specified
  if (!is.na(intervals[1])) {
    seq_remove <- c()
      for (i in 1:length(intervals)) {
        tmp_int <- intervals[i] %>% str_split("-") %>% unlist %>% as.numeric
        tmp_seq <- seq(tmp_int[1],tmp_int[2])
        tmp_seq <- tmp_seq[-c(1,length(tmp_seq))]
        seq_remove <- c(seq_remove, tmp_seq)
      }
    month.cut <- month.cut[!(month.cut %in% seq_remove)]
  }
  
  # copied copied from SUMMER::getBirths
  period.1yr <- year.cut[2] - year.cut[1] == 1
  surveyyear <- surveyyear - 1900
  year.cut <- year.cut - 1900
  variables <- union(variables, strata)
  if (!is.na(age.truncate)) {
    trunc <- intersect(which(dat[, age]%%12 == 0), which(dat[, 
                                                             age] >= age.truncate))
    if (length(trunc) > 0) {
      dat[trunc, age] <- dat[trunc, age] + 5
    } 
  } 
  
  datnew <- dat[, variables]
  dat[, alive] <- tolower(dat[, alive])
  datnew$dob <- dat[, dob] + cmc.adjust
  datnew$survey_year <- surveyyear
  
  # obsStart will be date of birth plus calendar adjustment (if any)
  datnew$obsStart <- dat[, dob] + cmc.adjust
  
  # date of death (dod) is date of birth + age + calendar adjustment (if any)
  datnew$dod <- dat[, dob] + dat[, age] + cmc.adjust
  
  # obsStop is dat of interview plus calendar adjustment (if any)
  datnew$obsStop <- dat[, date.interview] + cmc.adjust
  
  # for those who are dead, obsStop is actually date of death, so adjust this
  datnew$obsStop[(dat[, alive] == "no") | (dat[,alive] == "0")] <- datnew$dod[(dat[, alive] == "no") | (dat[, alive] == "0")]
  
  # create death indicator "died"
  datnew$died <- (dat[, alive] == "no") | (dat[,alive] == "0")
  
  # add child id as id.new
  datnew$id.new <- 1:nrow(datnew)
  
  # compute age at obsStop (either age of Death or age at censoring)
  datnew$age_at_censoring <- datnew$obsStop - datnew$obsStart
  
  # compute t0
  # for those who died: lower bound of interval
  # for those who did not die: max(age at censoring, 60 months)
  
  whichnotdied <- which(datnew$died == 0)
  # EDIT - fix this so that the interval ending with Inf doesn't have a closed right bracket
  datnew$age_interval <- cut(datnew$age_at_censoring, breaks = c(0,month.cut), include.lowest = TRUE, right = FALSE)
  
  lower_bound <- sapply(datnew$age_interval, function(x) {x %>% str_split(",") %>% unlist %>% nth(1) %>% 
      gsub("\\[|\\]", "", .) %>%
      gsub("\\(|\\)", "", .) %>%
      as.numeric})
  
  datnew$t0 <- lower_bound
  datnew[whichnotdied,]$t0 <- 0
  
  # compute t1
  # for those who died: upper bound of interval
  # for those who did not die: NA
  upper_bound <- sapply(datnew$age_interval, function(x) {x %>% str_split(",") %>% unlist %>% nth(2) %>% 
      gsub("\\[|\\]", "", .) %>%
      gsub("\\(|\\)", "", .) %>%
      as.numeric})
  
  datnew$t1 <- upper_bound
  datnew[whichnotdied,]$t1 <- Inf
  
  # make column for right_censored and interval_censored (these will be 1 - the other)
  datnew$right_censored <- !datnew$died
  datnew$interval_censored <- datnew$died
  
  # make column for year born (can base this on survey_year)
  # DHS recode says yyyy = int((CMC - 1) / 12) + 1900
  datnew$year_born <- floor((datnew$obsStart - 1)/12) + 1900
  datnew$month_born <- datnew$obsStart - ((datnew$year_born - 1900) * 12)
  
  # make column for year died
  datnew$year_died <- floor((datnew$dod - 1)/12) + 1900
  datnew$month_died <- datnew$dod - ((datnew$year_died - 1900) * 12)
  
  # fix up survey year
  datnew$survey_year <- datnew$survey_year + 1900
  
  # clean additional variable names
  datnew$cluster <- datnew$v001
  datnew$household <- datnew$v002
  
  if (length(strata) == 1) {
    datnew$strata <- datnew[,strata]
  } else {
    new_strata <- datnew[,strata[1]]
    for (i in 2:length(strata)) {
      new_strata <- paste0(new_strata, datnew[,strata[i]])
    }
    datnew$strata <- new_strata
  }
  
  datnew$urban <- datnew$v025
  datnew$weights <- datnew$v005
  
  # add in exact_age as recorded in B6 of DHS
  datnew$exact_age <- dat[,exact_age]
  
  # grab relevant columns
  datnew <- datnew[,c("caseid", "cluster", "household", "strata", "weights", "urban",
                      "survey_year", "id.new", "died", "t0", "t1", "age_at_censoring", "age_interval",
                      "year_born", "month_born", "year_died", "month_died", "right_censored",
                      "interval_censored", "exact_age")] 
  
  return(datnew)
}


#' Function to expand output from get_births into format appropriate for pseudolikelihood
#' 
#' @param surv_df output from get_births
#' @param period_boundaries a numeric vector containing the years surrounding each time period of interest.
#' Example: if periods are 2000-2005 and 2005-2010, \code{period_boundaries} should be \code{c(2000,2005,2010)}.
#' @param survey_date the (numeric) year of the DHS survey. If survey spand multiple years,
#' enter the latest one
#' 
#' @author Taylor Okonek
#' @noRd
#' @keywords internal
df_expand <- function(surv_df, 
                      period_boundaries,
                      survey_date) {
  
  # Get data in correct format
  p_boundaries <- period_boundaries
  p_boundaries <- p_boundaries %>% ISOdate(1, 1) %>% as.Date() %>% ymd()
  p <- 1:(length(p_boundaries) - 1)
  y_p <- p_boundaries[-length(p_boundaries)]
  l_p <- p_boundaries[2:length(p_boundaries)] - y_p
  survey_year <- survey_date 
  
  surv_df$individual <- 1:nrow(surv_df)
  
  df <- rcpp_expand_surv(surv_df, p, l_p, y_p)
  df$survey_year <- survey_year
  df$b_i <- df$b_i %>% ym()
  df$y_p <- df$y_p %>% ymd()
  df$a_pi <- round((df$y_p - df$b_i)/365 * 12) %>% as.numeric # age in months at start of time period
  
  df <- suppressMessages(left_join(df, surv_df[,c("individual","exact_age")]))
  
  surv_df_expanded <- df
  
  message("Data frame expanded. Getting indicator for individuals censored across time period boundaries")
  
  # convert l_p to months
  surv_df_expanded$l_p <- floor((surv_df_expanded$l_p %>% as.numeric())/30) 
  
  # Get A_i indicator for interval censored observations that are interval censored across
  # the boundary of a time period
  surv_df_expanded$A_i <- NA
  
  test <- rcpp_expand_surv_helper(surv_df_expanded, 
                                  nindiv = length(unique(surv_df_expanded$individual)),
                                  np = length(p))
  
  return(test)
  
}
