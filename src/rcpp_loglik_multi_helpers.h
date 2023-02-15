//sum.cpp
#include <Rcpp.h>
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
#include <boost/math/special_functions/gamma.hpp>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

NumericMatrix testDFtoNM1(DataFrame x) {
  int nRows=x.nrows();  
  NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
    y(_,i)=NumericVector(x[i]);
  }  
  return y;
}

// distributions:
// 0 = Exponential
// 1 = Weibull
// 2 = Piecewise Exponential
// [[Rcpp::export]]
double rcpp_hazard_integral(double lower_bound, double upper_bound, double log_shape, NumericVector log_scale_vec, int dist, NumericVector breakpoints) {
  NumericVector rate_param_vec = 1/exp(log_scale_vec);
  double shape_param = exp(log_shape);
  double ret_val = 0;
  int num_true = 0;
  LogicalVector U_p_internal(breakpoints.length());


  if (dist == 0 | dist == 1) {
    
    double rate_param = rate_param_vec[0];
    ret_val = pow(rate_param, shape_param) * (pow(upper_bound, shape_param) - pow(lower_bound, shape_param));

  } else if (dist == 2) {

    // are breakpoints relevant in this interval
    U_p_internal = (lower_bound < breakpoints) & (upper_bound > breakpoints);

    // get how many are true
    for (int i = 0; i < U_p_internal.length(); i++) {
      if (U_p_internal[i]) {
        num_true += 1;
      }
    }

    // get which values are true
    int which_are_true_counter = 0;
    NumericVector which_are_true(num_true);
    for (int i = 0; i < U_p_internal.length(); i++) {
      if (U_p_internal[i]) {
        which_are_true[which_are_true_counter] = breakpoints[i];
        which_are_true_counter += 1;
      }
    }

    // initialize new vector containing relevant breakpoints and lower_bound and upper_bound
    NumericVector new_bounds(2 + num_true);

    // fill in new_bounds, concatenate lower_bound, upper_bound, which_are_true
    new_bounds[0] = lower_bound;
    new_bounds[1] = upper_bound;
    if (which_are_true.length() > 0) {
      for (int i = 0; i < which_are_true.length(); i++) {
        new_bounds[2 + i] = which_are_true[i];
      }
    }

    // get unique bounds
    NumericVector new_bounds_unique = sort_unique(new_bounds);

    // get differences for each parameter
    NumericVector diffs(new_bounds_unique.length() - 1);
    for (int i = 0; i < diffs.length(); i++) {
      diffs[i] = new_bounds_unique[i + 1] - new_bounds_unique[i];
    }

    // get breakpoint intervals
    NumericVector breakpoints_intervals(breakpoints.length() + 2);
    NumericVector breakpoints_intervals_lower(breakpoints.length() + 1);
    NumericVector breakpoints_intervals_upper(breakpoints.length() + 1);

    breakpoints_intervals[0] = 0;
    breakpoints_intervals_lower[0] = 0;
    for (int i = 0; i < breakpoints.length(); i++) {
      breakpoints_intervals[i + 1] = breakpoints[i];
      breakpoints_intervals_lower[i + 1] = breakpoints[i];
      breakpoints_intervals_upper[i] = breakpoints[i];
    }
    breakpoints_intervals[breakpoints_intervals.length() - 1] = R_PosInf;
    breakpoints_intervals_upper[breakpoints_intervals_upper.length() - 1] = R_PosInf;

    // get relevant parameter ids
    LogicalVector relevant_pars = (lower_bound < breakpoints_intervals_upper) & (upper_bound > breakpoints_intervals_lower);

    // add up appropriate hazards
    int diff_counter = 0;
    for (int i = 0; i < rate_param_vec.length(); i++) {
      if (relevant_pars[i]) {
        ret_val += rate_param_vec[i] * diffs[diff_counter];
        diff_counter += 1;
      }
    }


  }
  return(ret_val);
}

// log hazard
// 0 = Exponential / Piecewise exponential
// 1 = Weibull
double rcpp_l_hazard(double x, double log_shape, NumericVector log_scale_vec, int dist, NumericVector breakpoints) {
  NumericVector rate_param_vec = 1/exp(log_scale_vec);
  double shape_param = exp(log_shape);
  double ret_val = 0;

  if (dist == 0) {

    ret_val = log(rate_param_vec[0]);

  } else if (dist == 1) {

    ret_val = log(rate_param_vec[0]) + log(shape_param) + (shape_param - 1) * log(x);

  } else {
    // get which parameter corresponds to the age when they died
    LogicalVector temp = x > breakpoints;
    double which_age_bin = sum(temp); // -1 because indexing
    
    ret_val = log(rate_param_vec[which_age_bin]); 

  }
  return(ret_val);
}

// log pdf
// 0 = Exponential
// 1 = Weibull
// 2 = Piecewise exponential - the way this is implemented here will work the same as for Exponential since it will identify the appropriate period
double rcpp_l_pdf(double x, double log_shape, double log_scale, int dist) {
  double rate_param = 1/exp(log_scale);
  double shape_param = exp(log_shape);
  double ret_val;
  if (dist == 0 | dist == 1 | dist == 2) {
    ret_val = log(rate_param * shape_param) + (shape_param - 1) * log(rate_param * x) - pow(rate_param * x, shape_param);
  }
  return(ret_val);
}