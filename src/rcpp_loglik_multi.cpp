//sum.cpp
#include <Rcpp.h>
#include "rcpp_loglik_multi_helpers.h"
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
#include <boost/math/special_functions/gamma.hpp>
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
//    2: exactly observed
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
// 2 = Piecewise Exponential
//
// breakpoints: vector of breakpoints for piecewise exponential distribution (not used for other distributions)
// par_period_id: integer vector of length j = 1, ..., J containing ids for which period log_scales[j] belongs to (values 1:P, possibly repeated)

// [[Rcpp::export]]
NumericVector rcpp_loglik_multi(DataFrame x_df, int num_periods, NumericVector log_shapes, NumericVector log_scales, int dist, NumericVector breakpoints, NumericVector par_period_id) {

  // Convert dataframe to matrix
  NumericMatrix x = testDFtoNM1(x_df);

  // Extract I_i, A_i, indicators for columns of x that refer to a_pi and l_p
  NumericVector I_i = x( _ , 0 );
  NumericVector A_i = x( _ , 1 );  

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
  double H_i, H_i0, H_i1, log_h_t;
  NumericVector a_pi(num_periods);
  NumericVector l_p(num_periods);
  NumericMatrix a_pi_mat, l_p_mat;
  int i;
  int j;
  int num_lived = 0;
  int which_vals_counter = 0;
  NumericVector j_vec(par_period_id.length());
  NumericVector which_vals(par_period_id.length() / num_periods);
  
  for (i = 0; i < x.nrow(); i++) {

    a_pi_mat = x(Range(i,i), Range(a_pi_cols[0], a_pi_cols[num_periods - 1]));
    l_p_mat = x(Range(i,i), Range(l_p_cols[0], l_p_cols[num_periods - 1]));
    
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
          which_vals_counter = 0;
          for (int k = 0; k < par_period_id.length(); k++) {
            if(par_period_id[k] == (j + 1)) {
              which_vals[which_vals_counter] = log_scales[k];
              which_vals_counter += 1;
            } 
          }
          H_i += rcpp_hazard_integral(std::max(a_pi[j], 0.0), std::min(x(i, 2), a_pi[j] + l_p[j]), log_shape_internal[j], which_vals, dist, breakpoints);
        }
      }
      
      NumericVector temp = 1/exp(which_vals);
      
      ret_vec[i] = -H_i;
      
      // if interval censored...
    } else if (I_i[i] == 1) {
      // identify which periods this person lived in for t_0

      U_p = (a_pi > -l_p) & (a_pi < x(i, 3));
      
      // accumulate hazard for t_0
      H_i0 = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          which_vals_counter = 0;
          for (int k = 0; k < par_period_id.length(); k++) {
            if(par_period_id[k] == (j + 1)) {
              which_vals[which_vals_counter] = log_scales[k];
              which_vals_counter += 1;
            } 
          }
          H_i0 += rcpp_hazard_integral(std::max(a_pi[j], 0.0), std::min(x(i, 3), a_pi[j] + l_p[j]), log_shape_internal[j], which_vals, dist, breakpoints);
        }
      }
      
      // identify which periods this person lived in for t_1
      
      U_p = (a_pi > -l_p) & (a_pi < x(i, 4));
      
      // accumulate hazard for t_1
      H_i1 = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          which_vals_counter = 0;
          for (int k = 0; k < par_period_id.length(); k++) {
            if(par_period_id[k] == (j + 1)) {
              which_vals[which_vals_counter] = log_scales[k];
              which_vals_counter += 1;
            } 
          }
          H_i1 += rcpp_hazard_integral(std::max(a_pi[j], 0.0), std::min(x(i, 4), a_pi[j] + l_p[j]), log_shape_internal[j], which_vals, dist, breakpoints);
        }
      }
      
      ret_vec[i] = log(exp(-H_i0) - exp(-H_i1));

      // if observed exactly...
    } else {

      // identify which periods this person lived in

      U_p = (a_pi > -l_p) & (a_pi < x(i, 3));

      int which_died = 0;
      H_i = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          which_died = j;
          which_vals_counter = 0;
          for (int k = 0; k < par_period_id.length(); k++) {
            if(par_period_id[k] == (j + 1)) {
              which_vals[which_vals_counter] = log_scales[k];
              which_vals_counter += 1;
            } 
          }
          H_i += rcpp_hazard_integral(std::max(a_pi[j], 0.0), std::min(x(i, 3), a_pi[j] + l_p[j]), log_shape_internal[j], which_vals, dist, breakpoints);
        }
      }

      // calculate log(h(t))
      which_vals_counter = 0;
      for (int k = 0; k < par_period_id.length(); k++) {
        if(par_period_id[k] == (which_died + 1)) {
          which_vals[which_vals_counter] = log_scales[k];
          which_vals_counter += 1;
        } 
      }

      log_h_t = rcpp_l_hazard(x(i, 3), log_shape_internal[which_died], which_vals, dist, breakpoints);

      // force it to have the same behavior as the interval censored observations, so that optim works
      // -inf + inf should not be NA, just let it be -inf
      if (R_PosInf == H_i) {
        ret_vec[i] = R_NegInf;
      } else {
        ret_vec[i] = -H_i + log_h_t;
      }

       
      
      
    }  

  }
  
  return(ret_vec);
}