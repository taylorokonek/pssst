//sum.cpp
#include <Rcpp.h>
#include "internal_helpers.h"
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame rcpp_expand_surv(DataFrame surv_df, NumericVector p, NumericVector l_p, CharacterVector y_p) {
  
  // Create individual_id vector for return dataframe
  NumericVector indiv(surv_df.nrow() * p.length());
  NumericVector p_long(surv_df.nrow() * p.length());
  NumericVector l_p_long(surv_df.nrow() * p.length());
  CharacterVector y_p_long(surv_df.nrow() * p.length());
  int which_ind = 1;
  int counter = 1;
  for (int i = 0; i < indiv.length(); i++) {
    indiv[i] = which_ind;
    counter += 1;
    if (counter == (p.length() + 1)) {
      which_ind += 1;
      counter = 1;
    }
  }

  // Initalize household, cluster, weights, strata, and survey_year columns
  NumericVector household(indiv.length());
  NumericVector cluster(indiv.length());
  CharacterVector strata(indiv.length());
  NumericVector weights(indiv.length());
  NumericVector I_i(indiv.length());
  CharacterVector b_i(indiv.length());
  NumericVector t_i(indiv.length());
  NumericVector t_0i(indiv.length());
  NumericVector t_1i(indiv.length());


  // Set up return data frame
  DataFrame ret_df = DataFrame::create(Named("individual") = clone(indiv),
    Named("household") = household,
    Named("cluster") = cluster,
    Named("strata") = strata,
    Named("weights") = weights,
    Named("p") = p_long,
    Named("y_p") = y_p_long,
    Named("l_p") = l_p_long,
    Named("b_i") = b_i,
    Named("t_i") = t_i,
    Named("t_0i") = t_0i,
    Named("t_1i") = t_1i,
    Named("I_i") = I_i);
  
  // Extract vectors from surv_df
  NumericVector surv_df_household = surv_df["household"];
  NumericVector surv_df_cluster = surv_df["cluster"];
  CharacterVector surv_df_strata = surv_df["strata"];
  NumericVector surv_df_weights = surv_df["weights"];
  IntegerVector surv_df_year_born = surv_df["year_born"];
  IntegerVector surv_df_month_born = surv_df["month_born"];
  LogicalVector surv_df_died = surv_df["died"];
  NumericVector surv_df_age_at_censoring = surv_df["age_at_censoring"];
  NumericVector surv_df_t0 = surv_df["t0"];
  NumericVector surv_df_t1 = surv_df["t1"];

  // Loop through individuals, and update vectors accordingly
  LogicalVector rel_rows(indiv.length());
  IntegerVector row_ids = seq_len(indiv.length()) - 1;
  NumericVector which_rows(p.length());
  for (int i = 0; i < surv_df.nrow(); i++) {
    rel_rows = indiv == (i + 1);
    which_rows = row_ids[rel_rows];

    for (int j = 0; j < which_rows.length(); j++) {
      household(which_rows[j]) = surv_df_household[i];
      cluster(which_rows[j]) = surv_df_cluster[i];
      strata(which_rows[j]) = surv_df_strata[i];
      weights(which_rows[j]) = surv_df_weights[i];
      p_long(which_rows[j]) = p[j];
      y_p_long(which_rows[j]) = y_p[j];
      l_p_long(which_rows[j]) = l_p[j];
      // Rcout << "here\n";
      b_i(which_rows[j]) = Rcpp::as<std::string> (paste3( Rcpp::as<std::string>( paste3(std::to_string((surv_df_year_born[i])), "-") ) , std::to_string((surv_df_month_born[i]))) );
      // Rcout << "here 2\n";

      if (surv_df_died[i]) {
        I_i(which_rows[j]) = 1;
        t_i(which_rows[j]) = R_NaN;
        t_0i(which_rows[j]) = surv_df_t0[i];
        t_1i(which_rows[j]) = surv_df_t1[i];
      } else {
        I_i(which_rows[j]) = 0;
        t_i(which_rows[j]) = surv_df_age_at_censoring[i];
        t_0i(which_rows[j]) = R_NaN;
        t_1i(which_rows[j]) = R_NaN;
      }
      
    }

  }


  return(ret_df);
}