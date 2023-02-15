//sum.cpp
#include <Rcpp.h>
#include <algorithm> 
#include <string.h>
#include <Rmath.h>
// #include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;

// rcpp_turnbull()
// 
// Assumes x_df has columns in this order: I_i, A_i, t_i, t_0i, t_1i, a_pi_1, ..., a_pi_p, l_p_1, ..., l_p_p
// Note: Currently only works with interval and right-censored data, does not deal with left-censored or exactly observed deaths
//
// niter: number of iterations to run algorithm
// t0: vector of lefthand side of intervals for censoring
// t1: vector of righthand side of intervals for censoring
// lefttrunc: vector of times at left truncation (0 if not left truncated)
// righttrunc: vector of times at right truncation (Inf if not right truncated)
// weights: vector of survey weights
// set_lower: vector of set lower bounds for intervals. Useful if there are set intervals you want to estimate that are more granular than your data
// set_upper: vector of set lower bounds for intervals. Useful if there are set intervals you want to estimate that are more granular than your data
//

// [[Rcpp::export]]
List rcpp_turnbull(int niter, NumericVector t0, NumericVector t1, NumericVector lefttrunc, NumericVector righttrunc, NumericVector weights, NumericVector set_lower, NumericVector set_upper) {

	// get intervals

	// first combine c(0, t0, t2, Inf)
	NumericVector temp(2 + t0.length() + t1.length());
	temp(0) = 0;
	temp(1) = R_PosInf;
	for (int i = 2; i < t0.length(); i++) {
		temp(i) = t0(i - 2);
		temp(i + t1.length()) = t1(i - 2);
	}
	NumericVector temp2 = sort_unique(temp);

	// declare q_j and p_j for lower and upper bounds of intervals
	NumericVector q_j, p_j;

	if (Rcpp::traits::is_nan<REALSXP>(set_lower[0])) {
		q_j  = NumericVector(temp2.length() - 1);
		p_j  = NumericVector(temp2.length() - 1);

		for (int i = 0; i < q_j.length(); i++) {
			q_j(i) = temp2(i);
			p_j(i) = temp2(i + 1);
		}
	} else {
		q_j  = NumericVector(set_lower.length());
		p_j  = NumericVector(set_upper.length());

		for (int i = 0; i < q_j.length(); i++) {
			q_j(i) = set_lower(i);
			p_j(i) = set_upper(i);
		}
	}


	// if anyone is exactly observed, make a closed interval [x,x]
	int num_closed = 0;
	LogicalVector same = (t0 == t1);
	for (int i = 0; i < same.length(); i++) {
		if (same(i)) {
			num_closed += 1;
		}
	}

	NumericVector which_closed(num_closed);
	if (num_closed > 0) {
		// figure out which values are closed intervals
		int j = 0;
		for (int i = 0; i < same.length(); i++) {
			if (same(i)) {
				which_closed(j) = t0(i);
				j += 1;
			}
		}
	}
	
	// get unique, exactly observed observations
	NumericVector which_closed_unique = sort_unique(which_closed);

	// update q_j, p_j to include which_closed
	NumericVector q_j_final, p_j_final;
	if (Rcpp::traits::is_nan<REALSXP>(set_lower[0])) {
	  q_j_final = NumericVector(q_j.length()  + which_closed_unique.length());
	  p_j_final = NumericVector(p_j.length()  + which_closed_unique.length());
	} else {
	  q_j_final = NumericVector(q_j.length());
	  p_j_final = NumericVector(p_j.length());
	}
	

	for (int i = 0; i < q_j.length(); i++) {
		q_j_final(i) = q_j(i);
		p_j_final(i) = p_j(i);
	}
	if (Rcpp::traits::is_nan<REALSXP>(set_lower[0])) {
	  for (int i = 0; i < which_closed_unique.length(); i++) {
	    q_j_final(i + q_j.length()) = which_closed_unique(i);
	    p_j_final(i + p_j.length()) = which_closed_unique(i);
	  }
	}

	q_j_final = q_j_final.sort();
	p_j_final = p_j_final.sort();

	// make lists of indicators for each individual
	// alpha_ij: indicator that the interval is a subset of their censoring set
	// beta_ij: indicator that the interval is a subset of their truncating set
	NumericMatrix alpha_ij(t0.length(), q_j_final.length());
	NumericMatrix beta_ij(t0.length(), q_j_final.length());
	double lower, upper, left, right;
	for (int i = 0; i < alpha_ij.nrow(); i++) {
		// get this person's censoring interval
		lower = t0(i);
		upper = t1(i);

		// get this person's truncation interval
		left = lefttrunc(i);
		right = righttrunc(i);

		alpha_ij( i , _ ) = (lower <= q_j_final) & (upper >= p_j_final);
		beta_ij( i , _ ) = (left <= q_j_final) & (right >= p_j_final);
	}

	// initialize s_j
	// mu_ij: probability i'th observation lies in [q_j, p_j]
	// nu_ij: expectation of number in the group corresponding to the i'th observation which have values in [q_j, p_j]
	NumericVector s_j(p_j_final.length());
	for (int i = 0; i < s_j.length(); i++) {
		s_j(i) = 1.0 / p_j_final.length();
	}

	NumericMatrix mu_ij(t0.length(), q_j_final.length());
	NumericMatrix nu_ij(t0.length(), q_j_final.length());

	// set up matrix to save s_j at each iteration
	NumericMatrix pi_j_mat(s_j.length(), niter);
	double M;
	NumericVector row_sums(mu_ij.nrow());
	NumericVector tempvec1(q_j_final.length());
	NumericVector tempvec2(q_j_final.length());
	NumericVector pi_j(s_j.length());

	for (int j = 0; j < niter; j++) {
		for (int i = 0; i < t0.length(); i++) {
			// get this person's censoring interval
			tempvec1 = alpha_ij(i,_) * s_j;
			tempvec2 = beta_ij(i,_) * s_j;
			mu_ij(i, _) = alpha_ij(i,_) * s_j / sum(tempvec1);
			nu_ij(i, _) = (1.0 - beta_ij(i,_)) * s_j / sum(tempvec2);
			tempvec1 = mu_ij(i, _);
			tempvec2 = nu_ij(i, _);
			row_sums(i) = sum(tempvec1 + tempvec2);
		}

		// calculate M and pi_j
		M = sum(weights * row_sums);

		// calculate pi_j
		for (int i = 0; i < s_j.length(); i++) {
			pi_j(i) = sum((mu_ij(_, i) + nu_ij(_, i)) * weights) / M;
		}

		pi_j_mat(_, j) = pi_j;

		// set s_j == pi_j
		s_j = pi_j;
	}

	// set up return list
	List ret_lst(3);

	ret_lst(0) = q_j_final;
	ret_lst(1) = p_j_final;
	ret_lst(2) = pi_j_mat;


	return(ret_lst);
}









