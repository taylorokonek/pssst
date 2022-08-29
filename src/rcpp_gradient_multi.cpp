//sum.cpp
#include <Rcpp.h>
#include "rcpp_gradient_multi_helpers.h"
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// Only works for interval- and right-censored data, no exactly observed times or left-censoring
// distributions:
// 0 = Exponential
// 1 = Weibull

// [[Rcpp::export]]
NumericVector rcpp_gradient_multi(DataFrame x_df, NumericVector log_shapes, NumericVector log_scales, int dist) {
    // Convert dataframe to matrix
  NumericMatrix x = testDFtoNM1_2(x_df);
  
  // Extract I_i, A_i, indicators for columns of x that refer to a_pi and l_p
  NumericVector I_i = x( _ , 0 );
  NumericVector A_i = x( _ , 1 );

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
  bool single_shape = FALSE;
  if (log_shapes.length() == 1) {
    single_shape = TRUE;
    for (int i = 0; i < num_periods; i++) {
      log_shape_internal(i) = log_shapes(0);
    }
  } else {
    log_shape_internal = log_shapes;
  }
  
  // initalize objects
  NumericMatrix ret_mat(x.nrow(), (num_periods + 1)*single_shape + (num_periods * 2)*(1 - single_shape));
  NumericVector grad_t0((num_periods + 1)*single_shape + (num_periods * 2)*(1 - single_shape));
  NumericVector grad_t1((num_periods + 1)*single_shape + (num_periods * 2)*(1 - single_shape));
  LogicalVector U_p;
  double H_i0, H_i1;
  NumericVector a_pi(num_periods);
  NumericVector l_p(num_periods);
  
  NumericMatrix a_pi_mat, l_p_mat;
  int i;
  int j;
  
  for (i = 0; i < x.nrow(); i++) {

    a_pi_mat = x(Range(i,i), Range(a_pi_cols[0], a_pi_cols[num_periods - 1]));
    l_p_mat = x(Range(i,i), Range(l_p_cols[0], l_p_cols[num_periods - 1]));
    
    a_pi = a_pi_mat( 0, _ );
    l_p = l_p_mat( 0, _ );
    
    // if right-censored...
    if (I_i[i] == 0) {
      // identify which periods this person lived in

      U_p = (a_pi > -l_p) & (a_pi < x(i, 2));

      // if only one shape parameter, gradient is a bit more complex
      if (single_shape) {
        // evaluate gradient with respect to log_shape
        ret_mat(i, 0) = - hazard_gradient_log_shape(x(i, 2), a_pi[U_p], l_p[U_p], log_shape_internal[0], log_scales[U_p], dist);

      // evaluate gradient with respect to log_scale
        for (int j = 0; j < log_scales.length(); j++) {
          if (U_p[j]) {
            ret_mat(i, j+1) = - hazard_gradient_log_scale(x(i, 2), a_pi[j], l_p[j], log_shape_internal[0], log_scales[j], dist);
          } else {
            ret_mat(i, j+1) = 0;
          }
        }

      // if multiple shape parameters, gradient is more simple
      } else {

        // evaluate gradient with respect to log_shape
        for (int j = 0; j < log_shape_internal.length(); j++) {
          if (U_p[j]) {
            ret_mat(i, j) = - hazard_gradient_log_shape_multi(x(i, 2), a_pi[j], l_p[j], log_shape_internal[j], log_scales[j], dist);
          } else {
            ret_mat(i, j) = 0;
          }
        }

        // evaluate gradient with respect to log_scale
        for (int j = 0; j < log_scales.length(); j++) {
          if (U_p[j]) {
            ret_mat(i, j + log_shape_internal.length()) = - hazard_gradient_log_scale(x(i, 2), a_pi[j], l_p[j], log_shape_internal[j], log_scales[j], dist);
          } else {
            ret_mat(i, j + log_shape_internal.length()) = 0;
          }
        }

      }
      
      
      // if interval censored...
    } else {
      // identify which periods this person lived in for t_0

      U_p = (a_pi > -l_p) & (a_pi < x(i, 3));
      
      // accumulate hazard for t_0
      H_i0 = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          H_i0 += rcpp_hazard_integral_2(std::max(a_pi[j], 0.0), std::min(x(i, 3), a_pi[j] + l_p[j]), log_shape_internal[j], log_scales[j], dist);
        }
      }

      if (single_shape) {
        // evaluate gradient with respect to log_shape at t_0
        grad_t0(0) = hazard_gradient_log_shape(x(i, 3), a_pi[U_p], l_p[U_p], log_shape_internal[0], log_scales[U_p], dist);

      // evaluate gradient with respect to log_scale at t_0
        for (int j = 0; j < log_scales.length(); j++) {
          if (U_p[j]) {
            grad_t0(j+1) = hazard_gradient_log_scale(x(i, 3), a_pi[j], l_p[j], log_shape_internal[0], log_scales[j], dist);
          } else {
            grad_t0(j+1) = 0;
          }
        }
      } else {

         // evaluate gradient with respect to log_shape at t_0
        for (int j = 0; j < log_shape_internal.length(); j++) {
          if (U_p[j]) {
            grad_t0(j) = hazard_gradient_log_shape_multi(x(i, 3), a_pi[j], l_p[j], log_shape_internal[j], log_scales[j], dist);
          } else {
            grad_t0(j) = 0;
          }
        }

      // evaluate gradient with respect to log_scale at t_0
        for (int j = 0; j < log_scales.length(); j++) {
          if (U_p[j]) {
            grad_t0(j+log_shape_internal.length()) = hazard_gradient_log_scale(x(i, 3), a_pi[j], l_p[j], log_shape_internal[j], log_scales[j], dist);
          } else {
            grad_t0(j+log_shape_internal.length()) = 0;
          }
        }

      }

      
      // identify which periods this person lived in for t_1
      
      U_p = (a_pi > -l_p) & (a_pi < x(i, 4));
      
      // accumulate hazard for t_1
      H_i1 = 0;
      for (j = 0; j < U_p.length(); j++) {
        if (U_p[j]) {
          H_i1 += rcpp_hazard_integral_2(std::max(a_pi[j], 0.0), std::min(x(i, 4), a_pi[j] + l_p[j]), log_shape_internal[j], log_scales[j], dist);
        }
      }

      
      if (single_shape) {
        // evaluate gradient with respect to log_shape at t_0
        grad_t1(0) = hazard_gradient_log_shape(x(i, 4), a_pi[U_p], l_p[U_p], log_shape_internal[0], log_scales[U_p], dist);

      // evaluate gradient with respect to log_scale at t_0
        for (int j = 0; j < log_scales.length(); j++) {
          if (U_p[j]) {
            grad_t1(j+1) = hazard_gradient_log_scale(x(i, 4), a_pi[j], l_p[j], log_shape_internal[0], log_scales[j], dist);
          } else {
            grad_t1(j+1) = 0;
          }
        }
      } else {

         // evaluate gradient with respect to log_shape at t_0
        for (int j = 0; j < log_shape_internal.length(); j++) {
          if (U_p[j]) {
            grad_t1(j) = hazard_gradient_log_shape_multi(x(i, 4), a_pi[j], l_p[j], log_shape_internal[j], log_scales[j], dist);
          } else {
            grad_t1(j) = 0;
          }
        }

      // evaluate gradient with respect to log_scale at t_0
        for (int j = 0; j < log_scales.length(); j++) {
          if (U_p[j]) {
            grad_t1(j+log_shape_internal.length()) = hazard_gradient_log_scale(x(i, 4), a_pi[j], l_p[j], log_shape_internal[j], log_scales[j], dist);
          } else {
            grad_t1(j+log_shape_internal.length()) = 0;
          }
        }

      }

      ret_mat(i, _) = (- exp(-H_i0) * grad_t0 + exp(-H_i1) * grad_t1) / (exp(-H_i0) - exp(-H_i1));
    }
  }
  
  return(ret_mat);
}