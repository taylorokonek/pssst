//sum.cpp
#include <Rcpp.h>
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

CharacterVector paste3(CharacterVector lhs, CharacterVector rhs)
  {
    using proxy_t = internal::string_proxy<STRSXP>;

    std::vector<std::string> res(lhs.begin(), lhs.end());
    std::transform(res.begin(), res.end(), rhs.begin(), res.begin(),
     [&](const std::string& x, const proxy_t& y) {
       return x + y;
     }
     );

    return wrap(res);
  }

double hazard_gradient_log_scale(double x, double a_pi, double l_p, double log_shape, double log_scale) {
  double beta_p = 1/exp(log_scale);
  double k = exp(log_shape);
    // beta_p = 1/scale_p
    double ret_val = k * pow(beta_p, k - 1) * (pow(std::min(x, a_pi + l_p),k) - pow(std::max(a_pi, 0.0),k)) * (-beta_p); // times -beta_p at the end because of chain rule
    return(ret_val);
  }

  double hazard_gradient_log_shape(double x, NumericVector a_pi, NumericVector l_p, double log_shape, NumericVector log_scales) {
    NumericVector beta_p = 1/exp(log_scales);
    double k = exp(log_shape);
    double ret_val = 0;

    for (int i = 0; i < a_pi.length(); i++) {
      ret_val += pow(beta_p[i] * std::min(x, a_pi[i] + l_p[i]), k) * log(beta_p[i] * std::min(x + 0.000001, a_pi[i] + l_p[i])) - pow(beta_p[i] * std::max(a_pi[i], 0.0), k) * log(beta_p[i] * std::max(a_pi[i], 0.000001));
    }
    return(ret_val * k); // times k here because of chain rule
  }

double hazard_gradient_log_shape_multi(double x, double a_pi, double l_p, double log_shape, double log_scales) {
    double beta_p = 1/exp(log_scales);
    double k = exp(log_shape);
    double ret_val;

    ret_val = pow(beta_p * std::min(x, a_pi + l_p), k) * log(beta_p * std::min(x + 0.000001, a_pi + l_p)) - pow(beta_p * std::max(a_pi, 0.0), k) * log(beta_p * std::max(a_pi, 0.000001));
    
    return(ret_val * k); // times k here because of chain rule
  }

  double rcpp_p_weibull(double x, double log_shape, double log_scale) {
    double ret_val;
    ret_val = 1 - exp(-pow((x/exp(log_scale)), exp(log_shape)));
    return(ret_val);
  }

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


  double rcpp_hazard_weibull(double x, double log_shape, double log_scale) {
    double rate_param = 1/exp(log_scale);
    double shape_param = exp(log_scale);
    return(rate_param * shape_param * pow(rate_param * x, shape_param - 1));
  }

  double rcpp_log_hazard_weibull(double x, double log_shape, double log_scale) {
    double rate_param = 1/exp(log_scale);
    double shape_param = exp(log_scale);

    return(log(shape_param - 1) + shape_param * log(rate_param) + (shape_param - 1) * log(x));

    // return(log(rate_param * shape_param) + (shape_param - 1) * log(rate_param * x));
  }

// set seed
  void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
  }