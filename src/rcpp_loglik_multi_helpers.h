//sum.cpp
#include <Rcpp.h>
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
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
double rcpp_hazard_integral(double lower_bound, double upper_bound, double log_shape, double log_scale, int dist) {
  double rate_param = 1/exp(log_scale);
  double shape_param = exp(log_shape);
  double ret_val;
  if (dist == 0 | dist == 1) {
    ret_val = pow(rate_param, shape_param) * (pow(upper_bound, shape_param) - pow(lower_bound, shape_param));
  }
  return(ret_val);
}

// log pdf
// 0 = Exponential
// 1 = Weibull
double rcpp_l_pdf(double x, double log_shape, double log_scale, int dist) {
  double rate_param = 1/exp(log_scale);
  double shape_param = exp(log_shape);
  double ret_val;
  if (dist == 0 | dist == 1) {
    ret_val = log(rate_param * shape_param) + (shape_param - 1) * log(rate_param * x) - pow(rate_param * x, shape_param);
  }
  return(ret_val);
}