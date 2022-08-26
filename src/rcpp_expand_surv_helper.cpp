//sum.cpp
#include <Rcpp.h>
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame rcpp_expand_surv_helper(DataFrame surv_df, double nindiv, double np) {
  // Loop through individuals, and update vectors accordingly
  LogicalVector rel_rows(surv_df.nrow());
  IntegerVector row_ids = seq_len(surv_df.nrow()) - 1;
  IntegerVector which_rows(np);

  // Extract columns from surv_df
  NumericVector surv_df_individual = surv_df["individual"];
  IntegerVector surv_df_I_i = surv_df["I_i"];
  IntegerVector surv_df_A_i = surv_df["A_i"];
  NumericVector surv_df_t0 = surv_df["t_0i"];
  NumericVector surv_df_t1 = surv_df["t_1i"];
  NumericVector surv_df_a_pi = surv_df["a_pi"];

  NumericVector temp_t0(np);
  NumericVector temp_a_pi(np);
  NumericVector temp_t1(np);
  LogicalVector temp(np);
  int k;

  for (int i = 0; i < nindiv; i ++) {
    rel_rows = surv_df_individual == (i + 1);
    which_rows = row_ids[rel_rows];

    // if they are right-censored, they can't have A_i == 1
    k = surv_df_I_i[which_rows[0]];

    if (k == 0) {
      for (int j = 0; j < np; j++) {
        surv_df_A_i(which_rows[j]) = 0;
      }

    // if they are interval-censored,
    } else {

      for (int j = 0; j < np; j++) {
        temp_t0(j) = surv_df_t0(which_rows[j]);
        temp_a_pi(j) = surv_df_a_pi(which_rows[j]);
        temp_t1(j) = surv_df_t1(which_rows[j]);
      }

      temp = (temp_t0 < temp_a_pi) & (temp_t1 > temp_a_pi);

      if (is_true(any(temp == TRUE))) {
        for (int j = 0; j < np; j++) {
          surv_df_A_i(which_rows[j]) = 1;
        }
      } else {
        for (int j = 0; j < np; j++) {
          surv_df_A_i(which_rows[j]) = 0;
        }
      }

    }
  }

  surv_df["A_i"] = surv_df_A_i;

  return(surv_df);

}