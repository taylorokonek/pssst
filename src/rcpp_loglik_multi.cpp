//sum.cpp
#include <Rcpp.h>
#include "rcpp_loglik_multi_helpers.h"
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// rcpp_loglik_multi()
// 
// Assumes x_df has columns in this order: I_i, A_i, t_i, t_0i, t_1i, a_pi_1, ..., a_pi_p, l_p_1, ..., l_p_p
// Note: Currently only works with interval and right-censored data, does not deal with left-censored or exactly observed deaths
//
// I_i: Indicator for which type of censoring we're dealing with:
//    0: right-censored at time t_i
//    1: interval-censored 
// A_i: Indicator for if they're interval censored across the boundary of a time period (can only be 1 if I_i = 1)
// t_i: Time at event
// t_0i: Left end point for interval censoring
// t_1i: Right end point for interval censoring
// a_pi_p: Age at beginning of period p
// l_p_p: Length of time of interval p  
//
// distributions:
// 0 = Exponential
// 1 = Weibull

// [[Rcpp::export]]
NumericVector rcpp_loglik_multi(DataFrame x_df, NumericVector log_shapes, NumericVector log_scales, int dist) {

  // Convert dataframe to matrix
  NumericMatrix x = testDFtoNM1(x_df);
  
  // Extract I_i, A_i, indicators for columns of x that refer to a_pi and l_p
  NumericVector I_i = x( _ , 0 );
  NumericVector A_i = x( _ , 1 );

  // Extract numer of periods (length of log_scales)
  double num_periods = log_scales.length();

  NumericVector a_pi_cols(num_periods);
  NumericVector l_p_cols(num_periods);
  for (int i = 0; i < num_periods; i++) {
    a_pi_cols[i] = 5 + i;
    l_p_cols[i] = 5 + num_periods + i;
  }
  
  // obtain parameters needed
  NumericVector log_shape_internal(num_periods);
  log_shape_internal(0) = log_shapes(0);
  if (log_shapes.length() == 1) {
    for (int i = 0; i < num_periods; i++) {
      log_shape_internal(i) = log_shapes(0);
    }
  } else {
    log_shape_internal = log_shapes;
  }
  
  // initalize objects
  NumericVector ret_vec(x.nrow());
  LogicalVector U_p;
  double H_i, H_i0, H_i1;
  NumericVector a_pi(num_periods);
  NumericVector l_p(num_periods);
  NumericMatrix a_pi_mat, l_p_mat;
  int i;
  int j;
  
  for (i = 0; i < x.nrow(); i++) {

    // Rcout << "Here 1 \n";

    a_pi_mat = x(Range(i,i), Range(a_pi_cols[0], a_pi_cols[num_periods - 1]));
    l_p_mat = x(Range(i,i), Range(l_p_cols[0], l_p_cols[num_periods - 1]));
    
    // Rcout << "Here 2 \n";
    
    a_pi = a_pi_mat( 0, _ );
    l_p = l_p_mat( 0, _ );
    
    // if right-censored...
    if (I_i[i] == 0) {
      // identify which periods this person lived in

      U_p = (a_pi > -l_p) & (a_pi < x(i, 2));
      
      // accumulate hazard
      H_i = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          H_i += rcpp_hazard_integral(std::max(a_pi[j], 0.0), std::min(x(i, 2), a_pi[j] + l_p[j]), log_shape_internal[j], log_scales[j], dist);
        }
      }
      
      ret_vec[i] = -H_i;
      
      // if interval censored...
    } else  {
      // identify which periods this person lived in for t_0

      U_p = (a_pi > -l_p) & (a_pi < x(i, 3));
      
      // accumulate hazard for t_0
      H_i0 = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          H_i0 += rcpp_hazard_integral(std::max(a_pi[j], 0.0), std::min(x(i, 3), a_pi[j] + l_p[j]), log_shape_internal[j], log_scales[j], dist);
        }
      }
      
      // identify which periods this person lived in for t_1
      
      U_p = (a_pi > -l_p) & (a_pi < x(i, 4));
      
      // accumulate hazard for t_1
      H_i1 = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          H_i1 += rcpp_hazard_integral(std::max(a_pi[j], 0.0), std::min(x(i, 4), a_pi[j] + l_p[j]), log_shape_internal[j], log_scales[j], dist);
        }
      }
      
      ret_vec[i] = log(exp(-H_i0) - exp(-H_i1));

    }  

  }
  
  return(ret_vec);
}