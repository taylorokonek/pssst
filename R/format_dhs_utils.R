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
  # apply for intervals where everyone within the interval is interval-censored
  if (is.null(names(intervals))) {
    intervals_general <- intervals
    intervals_targeted <- NULL
  } else {
    intervals_general <- intervals[which(names(intervals)=="")]
    intervals_targeted <- intervals[which(names(intervals)!="")]
  }
  if (!all(is.na(intervals_general)) & length(intervals_general) != 0) {
    seq_remove <- c()
      for (i in 1:length(intervals_general)) {
        tmp_int <- intervals_general[i] %>% str_split("-") %>% unlist %>% as.numeric
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
  datnew$age_interval <- as.character(cut(datnew$age_at_censoring, breaks = c(0,month.cut), include.lowest = TRUE, right = FALSE))
  
  # extra interval censoring, for specific months
  if (length(intervals_targeted) > 0) {
    for (i in 1:length(intervals_targeted)) {
      censor_age <- as.numeric(names(intervals_targeted)[i])
      age_interval <- paste0("[", gsub("-", ",", unname(intervals_targeted)[i]), ")")
      datnew[!is.na(dat[,age]) & dat[,age] == censor_age,]$age_interval <- age_interval
    }
  }
  
  # use exact age interval for deaths recorded in years
  # deaths before 24 months are supposed to be recorded in months but sometimes they are recorded as 1 year
  # convert exact_column to numeric just for this one bit of code
  age_in_years <- which(as.numeric(dat[, exact_age]) >= 300 & as.numeric(dat[, exact_age]) < 350)
  if (length(age_in_years) > 0) {
    age_in_years_vals <- as.numeric(dat[age_in_years, exact_age])-300
    datnew[age_in_years,]$age_interval <- paste0("[", age_in_years_vals*12, ",", age_in_years_vals*12+12, ")")
  }
  
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


#' @title Add strata column to DHS birth recode module
#'
#' @description Add strata column for DHS birth recode manual, based on
#'  combining existing columns, e.g. a column for region and a column for
#'  rurality. This is based on the DHS .do file available at:
#'  https://github.com/DHSProgram/DHS-Analysis-Code/blob/main/Intro_DHSdata_Analysis/2_SurveyDesign/Survey_strata.do
#'  and is consistent with the guidance in Chapter 1 of DHS Guide to Statistics.
#'  
#' @param df data.frame with DHS dataset as downloaded.
#' 
#' @references
#' Assaf, S. (2021). Survey_strata.do. GitHub.
#'  https://github.com/DHSProgram/DHS-Analysis-Code/blob/main/Intro_DHSdata_Analysis/2_SurveyDesign/Survey_strata.do
#' 
#' Croft, Trevor N., Allen, Courtney K., Zachary, Blake W., et al. 2023.
#'  Guide to DHS Statistics. Rockville, Maryland, USA: ICF. https://dhsprogram.com/data/Guide-to-DHS-Statistics/index.cfm
#' 
#' @return data.frame with "strata" column added.
#'
#' @export
create_strata_dhs <- function(df) {
  
  # add missing columns as necessary, fill with NA
  cols <- c("v022", "v023", "v024", "v025", "sstratum", "sstrate", "sdepart",
            "scty", "secozone", "sdepto")
  df[setdiff(cols, names(df))] <- NA
  
  # create default strata and some other new variables
  df <- df %>%
    mutate(
      strata = as.character(v022), # default strata
      strata_v023 = paste0(v023, " - ", v025),
      strata_v024 = paste0(v024, " - ", v025), # v024 is region, v025 is rurality
      phase = v000,
      year = as.numeric(as.character(v007))
    )
  
  # list exception groups
  v023_phases <- 
    c("AL5", "AM4", "AZ5", "BD5", "BJ4", "BF3", "CM4", "KM6", "CI3", "EG2",
      "ET4", "GA3", "GH3", "GN3", "ID2", "ID3", "KK3", "KE2", "KE3", "KY3",
      "LB6", "MD3", "MD4", "MW2", "ML3", "MB4", "MZ3", "NM4", "NG3", "NG4",
      "PE5", "SN2", "TZ2", "TR3", "TR4", "TR6", "UG3", "UZ3", "ZM3", "ZW3",
      "ZW4")
  strata_v023_phases <- 
    c("BD3", "BO4", "BR2", "BR3", "BF4", "DR2", "GH4", "GU3", "GN4", "HN5",
      "JO5", "KE4", "MW4", "ML4", "ML5", "MA4", "MZ4", "NP4", "NP5", "NP6",
      "NC3", "NC4", "NI5", "PK5", "PK6", "PE2", "PE3", "PH2", "PH3", "PH4",
      "RW4", "RW6", "TR2")
  strata_v024_phases <- 
    c("BD4", "BJ3", "BJ6", "CF3", "TD4", "CO2", "CO3", "KM3", "CG5", "CD5",
      "DR3", "DR4", "DR5", "EG3", "EG4", "GH2", "GY4", "HT4", "HT5", "IA2",
      "IA3", "ID4", "ID5", "JO2", "JO3", "LS4", "LS5", "MD2", "MA2", "NM2",
      "NM5", "NI2", "NI3", "NG2", "PK2", "PY2", "PE4", "PH5", "RW2", "ST5",
      "SN4", "SZ5", "ZA3", "TZ3", "TG3", "TR5", "UG4", "UG5", "UG6", "UA5",
      "VNT", "YE2", "ZM2", "ZM4", "ZW", "ZW5", "SN5") 
  group_101_102_phases <- 
    c("BR", "BF2", "BU", "CO", "DR", "EC", "ES", "EG", "GH", "GU", "ID", "KE",
      "LB", "LK", "ML", "MX", "MA", "OS", "PE", "SD", "SN", "TH", "TT", "TG",
      "TN", "UG")
  
  # apply exceptions
  df <- df %>%
    mutate(
      strata = case_when(
        phase %in% v023_phases ~ as.character(v023),
        phase %in% strata_v023_phases ~ strata_v023,
        phase %in% strata_v024_phases ~ strata_v024,
        phase == "BD3" & year %in% c(99, 2000, 0) ~ as.character(v023),
        phase == "BO3" & year %in% c(93, 94) ~ strata_v023,
        phase %in% group_101_102_phases ~ paste0(v101, " - ", v102),
        phase == "KH4" ~ as.character(sstratum),
        phase == "KH5" & year %in% c(2005, 2006) ~ strata_v024,
        phase %in% c("CM2", "CM3") ~ as.character(sstrate),
        phase == "CO4" & year == 2000 ~ strata_v023,
        phase == "CO4" & year %in% c(2004, 2005) ~ strata_v024,
        phase == "HT3" ~ paste0(sdepart, " - ", v025),
        phase == "LB5" & year < 2008 ~ paste0(scty, " - ", v025, " - ", v024),
        phase == "NP3" ~ paste0(secozone, " - ", v024),
        phase == "SN2" & year %in% c(92, 93) ~ strata_v024,
        phase == "BO" & year == 89 ~ paste0(sdepto, " - ", v102),
        TRUE ~ strata
      )
    )
  
  return(df)
}

