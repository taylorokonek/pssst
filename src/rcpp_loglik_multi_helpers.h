//sum.cpp
#include <Rcpp.h>
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
#include <boost/math/special_functions/gamma.hpp>
#include <RcppNumerical.h>
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

// cribbed from https://gallery.rcpp.org/articles/vector-cumulative-sum/
NumericVector cumsum1(NumericVector x){
  // initialize an accumulator variable
  double acc = 0;
  
  // initialize the result vector
  NumericVector res(x.size());
  
  for(int i = 0; i < x.size(); i++){
    acc += x[i];
    res[i] = acc;
  }
  return res;
}

// input shape = Q, location = omega, scale = sigma
// [[Rcpp::export]]
double rcpp_F_gengamma(double x, double alpha, double beta, double gamma, bool lower_tail, bool give_log) {
  // double d = gamma * alpha;
  // double a = 1/beta;
  // double p = gamma;
  
  // // double numerator = boost::math::tgamma_lower(alpha, pow(x/a, gamma));
  // double numerator = boost::math::gamma_p(alpha, pow(x/a, gamma)); // use normalized/regularized lower incomplete gamma function
  double Q = alpha;
  double mu = beta;
  double sigma = gamma;
  double omega = -(mu - log(x)) / sigma;
  
  // double ret_val = boost::math::gamma_p(pow(Q, -2), pow(exp(-omega) * exp(-2 * sigma * log(Q) / Q) * x, Q / sigma));
  // double ret_val = R::pgamma(pow(exp(-omega) * exp(-2 * sigma * log(Q) / Q) * x, Q / sigma), pow(Q, -2), 1, 1, 0);
  double ret_val = R::pgamma(pow(Q, -2) * exp(Q * omega), pow(Q, -2), 1, lower_tail, give_log);
  
  // double denominator = tgamma(alpha);
  // double ret_val = numerator;
  
  return(ret_val);
}

// input shape = Q, location = omega, scale = sigma
// [[Rcpp::export]]
double rcpp_f_gengamma(double x, double alpha, double beta, double gamma, bool log_pdf) {
  double Q = alpha;
  double mu = beta;
  double sigma = gamma;
  
  Environment pkg = Environment::namespace_env("flexsurv");
  Function f = pkg["dgengamma"];
  SEXP temp = f(x, Named("mu") = mu, Named("Q") = Q, Named("sigma") = sigma, Named("log") = log_pdf);
  
  return Rcpp::as<double>(temp);
  
}

// [[Rcpp::export]]
double rcpp_F_gompertz(double x, double rate, double shape, bool lower_tail, bool give_log) {
  
  double a = shape;
  double b = rate;
  double ret_val;
  
  ret_val = -(b/a) * (exp(a * x) - 1);
  // ret_val = 1 - exp(-(b/a) * (exp(a * x) - 1));
  // ret_val = 1 - exp(-eta * (exp(b * x - 1)));
  
  if (lower_tail) {
    if (give_log) {
      ret_val = log1p(-exp(ret_val));
    } else {
      ret_val = -expm1(ret_val);
    }
  } else {
    if (give_log) {
      ret_val = ret_val;
    }  else {
      ret_val = exp(ret_val);
    }
  }
  
  return(ret_val);
}

// [[Rcpp::export]]
double rcpp_F_loglogistic(double x, double shape, double scale, bool lower_tail, bool give_log) {
  
  double alpha = scale;
  double beta = shape;
  double ret_val;
  
  ret_val = 1 - pow(1 + pow(x/alpha, beta), -1);
  
  if (lower_tail) {
    if (give_log) {
      ret_val = log(ret_val);
    } else {
      ret_val = ret_val;
    }
  } else {
    if (give_log) {
      ret_val = log(1 - ret_val);
    }  else {
      ret_val = 1 - ret_val;
    }
  }
  
  return(ret_val);
}

// [[Rcpp::export]]
double rcpp_f_loglogistic(double x, double shape, double scale, bool give_log) {
  double alpha = scale;
  double beta = shape;
  double ret_val;
  
  ret_val = (beta/alpha) * pow(x/alpha, beta - 1) / pow(1 + pow(x/alpha, beta), 2);
  
  
  if (give_log) {
    ret_val = log(ret_val);
  } else {
    ret_val = ret_val;
  }
  
  return(ret_val);
  
}

// [[Rcpp::export]]
double rcpp_F_dagum(double x, double shape1, double shape2, double scale, bool lower_tail, bool give_log) {
  double p = shape1;
  double a = shape2;
  double b = scale;
  double ret_val;
  ret_val = pow( 1 + pow(x/b, -1*a), -1*p);
  
  if (lower_tail) {
    if (give_log) {
      ret_val = log(ret_val);
    } else {
      ret_val = ret_val;
    }
  } else {
    if (give_log) {
      ret_val = log(1 - ret_val);
    }  else {
      ret_val = 1 - ret_val;
    }
  }
  return(ret_val);
}

// [[Rcpp::export]]
double rcpp_f_dagum(double x, double shape1, double shape2, double scale, bool give_log) {
  double p = shape1;
  double a = shape2;
  double b = scale;
  double ret_val;
  ret_val = ((a*p)/x)*(pow(x/b, a*p))/ pow((pow(x/b, a) + 1), p+1);
  if (give_log) {
    ret_val = log(ret_val);
  } else {
    ret_val = ret_val;
  }
  return(ret_val);
}



// input shape = Q, location = omega, scale = sigma
// [[Rcpp::export]]
double rcpp_f_gompertz(double x, double rate, double shape, bool give_log) {
  
  double a = shape;
  double b = rate;
  double ret_val;
  
  ret_val = log(b) + a * x - (b/a) * (exp(a * x) - 1);
  // ret_val = log(b) + log(eta) + eta + b * x - eta * exp(b * x);
  
  if (give_log) {
    ret_val = ret_val;
  } else {
    ret_val = exp(ret_val);
  }
  
  return(ret_val);
}

// Exponentially-truncated shifted-power family of hazards
class etsp_haz: public Numer::Func
{
private:
  double a;
  double b;
  double c;
  double p;
public:
  etsp_haz(double a_, double b_, double c_, double p_) : a(a_), b(b_), c(c_), p(p_) {}
  
  double operator()(const double& x) const
  {
    return a * pow(x + c, -p) * exp(-b * x);
  }
};

// distributions:
// 0 = Exponential
// 1 = Weibull
// 2 = Piecewise Exponential
// 3 = Generalized Gamma
// 4 = lognormal
// 5 = Gompertz
// 6 = ETSP: Exponentially-trunacted shifted-power
// 7 = Log-logistic
// 8 = Dagum
// [[Rcpp::export]]
double rcpp_hazard_integral(double lower_bound, double upper_bound, double log_shape, NumericVector log_scale_vec, int dist, NumericVector breakpoints, double etsp_c = 0.0, bool force_nonincreasing = false) {
  NumericVector scale_vec = exp(log_scale_vec);
  NumericVector rate_param_vec = 1/scale_vec;

  if (force_nonincreasing) {
    if (dist == 2) {
      rate_param_vec = cumsum1(rate_param_vec);
      std::reverse(rate_param_vec.begin(), rate_param_vec.end());
    }
  }

  double shape_param = exp(log_shape);

  if (force_nonincreasing) {
    if (dist == 7) {
      shape_param = exp(log_shape)/(1+exp(log_shape)); // logit instead of log transformation for log-logistic to constrain non-increasing hazard
    }
  }
  
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
    
  } else if (dist == 3) {
    
    // rate_param_vec = c(beta, alpha)
    // shape_param = gamma
    
    // ret_val = -log(1 - rcpp_F_gengamma(upper_bound, scale_vec[1], scale_vec[0], shape_param, 1, 0)) + log(1 - rcpp_F_gengamma(lower_bound, scale_vec[1], scale_vec[0], shape_param, 1, 0));
    ret_val = - rcpp_F_gengamma(upper_bound, scale_vec[1], scale_vec[0], shape_param, 0, 1) + rcpp_F_gengamma(lower_bound, scale_vec[1], scale_vec[0], shape_param, 0, 1);
    
  } else if (dist == 4) {
    
    double rate_param = rate_param_vec[0];
    ret_val = - R::plnorm(upper_bound, rate_param, shape_param, 0, 1) + R::plnorm(lower_bound, rate_param, shape_param, 0, 1);
    
  } else if (dist == 5) {
    double rate_param = rate_param_vec[0];
    ret_val = - rcpp_F_gompertz(upper_bound, rate_param, shape_param, 0, 1) + rcpp_F_gompertz(lower_bound, rate_param, shape_param, 0, 1);
    
  } else if (dist == 6) {
    
    // a = shape_param
    // b = scale_vec[0]
    // c = etsp_c
    // p = scale_vec[1]
    
    etsp_haz f(shape_param, scale_vec[0], etsp_c, scale_vec[1]);
    double err_est;
    int err_code;
    
    ret_val += Numer::integrate(f, lower_bound, upper_bound, err_est, err_code);
    
  } else if (dist == 7) {
    
    double scale_param = scale_vec[0];
    
    ret_val = - rcpp_F_loglogistic(upper_bound, shape_param, scale_param, 0, 1) + rcpp_F_loglogistic(lower_bound, shape_param, scale_param, 0, 1);
  }
  else if (dist == 8) {
    
    // the Dagum distribution has 3 parameters:
    // p and a are shape parameters, b is the scale parameter
    // p = shape1, a = shape2 
    // for compatibility with existing data structures, 
    // shape2 is set to be the first entry of the scale vector
    
    double shape1 = shape_param; 
    double shape2 = scale_vec[0];
    double scale = scale_vec[1];
    
    ret_val = - rcpp_F_dagum(upper_bound, shape1, shape2, scale, 0, 1) + rcpp_F_dagum(lower_bound, shape1, shape2, scale, 0, 1);
    
  }
  
  return(ret_val);
}

// log hazard
// 0 = Exponential
// 1 = Weibull
// 2 = Piecewise exponential
// 3 = Generalized gamma
// 4 = lognormal
// 5 = Gompertz
// 7 = Log-logistic
// 8 = Dagum
double rcpp_l_hazard(double x, double log_shape, NumericVector log_scale_vec, int dist, NumericVector breakpoints, double etsp_c = 0.0, bool force_nonincreasing = false) {
  NumericVector scale_vec = exp(log_scale_vec);
  NumericVector rate_param_vec = 1/scale_vec;

  if (force_nonincreasing) {
    if (dist == 2) {
      rate_param_vec = cumsum1(rate_param_vec);
      std::reverse(rate_param_vec.begin(), rate_param_vec.end());
    }
  }

  double shape_param = exp(log_shape);

  if (force_nonincreasing) {
    if (dist == 7) {
      shape_param = exp(log_shape)/(1+exp(log_shape)); // logit instead of log transformation for log-logistic to constrain non-increasing hazard
    }
  }

  double ret_val = 0;
  double numerator, denominator;
  
  if (dist == 0) {
    
    ret_val = log(rate_param_vec[0]);
    
  } else if (dist == 1) {
    
    ret_val = log(rate_param_vec[0]) + log(shape_param) + (shape_param - 1) * log(x);
    
  } else if (dist == 2) {
    // get which parameter corresponds to the age when they died
    LogicalVector temp = x > breakpoints;
    double which_age_bin = sum(temp); // -1 because indexing
    
    ret_val = log(rate_param_vec[which_age_bin]); 
    
  } else if (dist == 3) {
    // h(t) = f(t)/(1 - F(t))
    // rate_param_vec = c(beta, alpha)
    // shape_param = gamma
    numerator = rcpp_f_gengamma(x, scale_vec[1], scale_vec[0], shape_param, 1);
    denominator = rcpp_F_gengamma(x, scale_vec[1], scale_vec[0], shape_param, 0, 1);
    ret_val = numerator - denominator;
    
  } else if (dist == 4) {
    
    numerator = R::dlnorm(x, rate_param_vec[0], shape_param, 1);
    denominator = R::plnorm(x, rate_param_vec[0], shape_param, 0, 1);
    ret_val = numerator - denominator;
    
  } else if (dist == 5) {
    numerator = rcpp_f_gompertz(x, rate_param_vec[0], shape_param, 1);
    denominator = rcpp_F_gompertz(x, rate_param_vec[0], shape_param, 0, 1);
    ret_val = numerator - denominator;
    
  } else if (dist == 6) {
    double a = shape_param;
    double b = scale_vec[0];
    double c = etsp_c;
    double p = scale_vec[1];
    
    ret_val = log(a) - p * log(x + c) - b * x;
  } else if (dist == 7) {
    
    numerator = rcpp_f_loglogistic(x, shape_param, scale_vec[0], 1);
    denominator = rcpp_F_loglogistic(x, shape_param, scale_vec[0], 0, 1);
    ret_val = numerator - denominator;
    
  }
  else if (dist == 8) {
    
    numerator = rcpp_f_dagum(x, shape_param, scale_vec[0], scale_vec[1], 1);
    denominator = rcpp_F_dagum(x, shape_param, scale_vec[0], scale_vec[1], 0, 1);
    ret_val = numerator - denominator;
  }
  
  return(ret_val);
}

// log pdf
// 0 = Exponential
// 1 = Weibull
// 2 = Piecewise exponential - the way this is implemented here will work the same as for Exponential since it will identify the appropriate period
// 3 = Generalized gamma
// 4 = lognormal
// 5 = Gompertz
double rcpp_l_pdf(double x, double log_shape, double log_scale, int dist) {
  double rate_param = 1/exp(log_scale);
  double shape_param = exp(log_shape);
  double ret_val;
  if (dist == 0 | dist == 1 | dist == 2) {
    ret_val = log(rate_param * shape_param) + (shape_param - 1) * log(rate_param * x) - pow(rate_param * x, shape_param);
  }
  return(ret_val);
}