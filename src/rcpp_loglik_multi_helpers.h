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

  double rcpp_hazard_weibull_integral(double lower_bound, double upper_bound, double log_shape, double log_scale) {
    double rate_param = 1/exp(log_scale);
    double shape_param = exp(log_shape);
    return(pow(rate_param, shape_param) * (pow(upper_bound, shape_param) - pow(lower_bound, shape_param)));
  }