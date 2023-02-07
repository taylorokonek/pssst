//sum.cpp
#include <Rcpp.h>
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

NumericMatrix testDFtoNM1_2(DataFrame x) {
  int nRows=x.nrows();  
  NumericMatrix y(nRows,x.size());
  for (int i=0; i<x.size();i++) {
    y(_,i)=NumericVector(x[i]);
  }  
  return y;
}

double hazard_gradient_log_scale(double x, double a_pi, double l_p, double log_shape, double log_scale, int dist) {
  double beta_p = 1/exp(log_scale);
  double k = exp(log_shape);
  double ret_val;
  // beta_p = 1/scale_p
  if (dist == 0 | dist == 1) {
    ret_val = k * pow(beta_p, k - 1) * (pow(std::min(x, a_pi + l_p),k) - pow(std::max(a_pi, 0.0),k)) * (-beta_p); // times -beta_p at the end because of chain rule
  }
  return(ret_val);
}

double hazard_gradient_log_shape(double x, NumericVector a_pi, NumericVector l_p, double log_shape, NumericVector log_scales, int dist) {
  NumericVector beta_p = 1/exp(log_scales);
  double k = exp(log_shape);
  double ret_val = 0;
  
  if (dist == 0 | dist == 1) {
    for (int i = 0; i < a_pi.length(); i++) {
      ret_val += pow(beta_p[i] * std::min(x, a_pi[i] + l_p[i]), k) * log(beta_p[i] * std::min(x + 0.000001, a_pi[i] + l_p[i])) - pow(beta_p[i] * std::max(a_pi[i], 0.0), k) * log(beta_p[i] * std::max(a_pi[i], 0.000001));
    }
  }
  
  return(ret_val * k); // times k here because of chain rule
}

double hazard_gradient_log_shape_multi(double x, double a_pi, double l_p, double log_shape, double log_scales, int dist) {
  double beta_p = 1/exp(log_scales);
  double k = exp(log_shape);
  double ret_val;
  
  if (dist == 0 | dist == 1) {
    ret_val = pow(beta_p * std::min(x, a_pi + l_p), k) * log(beta_p * std::min(x + 0.000001, a_pi + l_p)) - pow(beta_p * std::max(a_pi, 0.0), k) * log(beta_p * std::max(a_pi, 0.000001));
    
  }
  
  return(ret_val * k); // times k here because of chain rule
}

double l_pdf_gradient_log_scale(double x, double log_shape, double log_scale, int dist) {
  double beta_p = 1/exp(log_scale);
  double k = exp(log_shape);
  double ret_val;
  // beta_p = 1/scale_p
  if (dist == 0 | dist == 1) {
    ret_val = -(k - k * pow(beta_p * x, k)); // times -beta_p because chain rule at end
  }
  return(ret_val); 
}

double l_pdf_gradient_log_shape(double x, double log_shape, double log_scale, int dist) {
  double beta_p = 1/exp(log_scale);
  double k = exp(log_shape);
  double ret_val;
  // beta_p = 1/scale_p
  if (dist == 0 | dist == 1) {
    ret_val = 1 - k * (pow(beta_p * x, k) - 1) * log(beta_p * x); // times k at the end because chain rule
  }
  return(ret_val); 
}

double rcpp_hazard_integral_2(double lower_bound, double upper_bound, double log_shape, double log_scale, int dist) {
  double rate_param = 1/exp(log_scale);
  double shape_param = exp(log_shape);
  double ret_val;
  if (dist == 0 | dist == 1) {
    ret_val = pow(rate_param, shape_param) * (pow(upper_bound, shape_param) - pow(lower_bound, shape_param));
  }
  return(ret_val);
}