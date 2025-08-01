% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surv_synthetic.R
\name{surv_synthetic}
\alias{surv_synthetic}
\title{Fit a parametric survival model to DHS data across multiple time periods.}
\usage{
surv_synthetic(
  df,
  individual = "individual",
  survey = TRUE,
  household = NULL,
  cluster = NULL,
  strata = NULL,
  weights = NULL,
  p = "p",
  a_pi = "a_pi",
  l_p = "l_p",
  I_i = "I_i",
  A_i = "A_i",
  t_i = "t_i",
  t_0i = "t_0i",
  t_1i = "t_1i",
  only_scale = FALSE,
  numerical_grad = FALSE,
  dist = "weibull",
  breakpoints = NA,
  init_vals = NA,
  etsp_c = 0,
  force_nonincreasing = FALSE
)
}
\arguments{
\item{df}{a dataframe containing the output from \code{format_dhs}, or optionally, dataframe
containing the following columns}

\item{individual}{column corresponding to individual ID in \code{df}}

\item{survey}{boolean for whether or not your data comes from a survey with a specified probability design.
If true, must specify household, cluster, strata, and weights parameters in the function, and a finite
population variance based on a pseudo-likelihood will be returned for your survey-weighted estimates.
If false, superpopulation estimates will be returned with the usual asymptotic variance estimator.
Defaults to \code{TRUE}.}

\item{household}{column corresponding to household ID in \code{df}}

\item{cluster}{column corresponding to cluster ID in \code{df}}

\item{strata}{column corresponding to strata ID in \code{df}}

\item{weights}{column corresponding to weights in \code{df}. If survey = FALSE, you may
still specify a weights column, but note that the superpopulation variance will be returned
as opposed to the finite population variance.}

\item{p}{column corresponding to period ID (must be integer-valued, numeric) in \code{df}}

\item{a_pi}{column corresponding to child's age at beginning of time period in \code{df}}

\item{l_p}{column corresponding to length of time period in \code{df}}

\item{I_i}{column corresponding to indicator for interval censoring in \code{df}. Should be 1 if
child is interval-censored, 0 if right-censored.}

\item{A_i}{column corresponding to indicator for if a child is interval-censored across the
boundary of a time period in \code{df}}

\item{t_i}{column corresponding to age at right-censoring, if right-censored, in \code{df}}

\item{t_0i}{column corresponding to lower bound of interval, if interval-censored, in \code{df}}

\item{t_1i}{column corresponding to upper bound of interval, if interval-censored, in \code{df}}

\item{only_scale}{boolean for varying only the scale parameter across time period. Defaults to
\code{FALSE}. This option is only available for dist = "weibull"}

\item{numerical_grad}{boolean for whether gradient should be calculated numerically or
analytically. Analytical gradient is faster, but only available for Weibull and Exponential distributions
at the moment.}

\item{dist}{distribution. Currently supports "weibull", "exponential",
"piecewise_exponential", "gengamma", "lognormal", "gompertz", "etsp" (exponentially-truncated shifted power family),
"loglogistic", "dagum". The loglogistic and piecewise exponential distributions can optionally be constrained
to have a monotonically non-increasing hazard, using \code{force_nonincreasing}. See Details for more.}

\item{breakpoints}{if distribution is "piecewise_exponential", the breakpoints (in months) where
the distribution should be divided. For example, to have three distinct age groups [0,1), [1,12), [12,\eqn{\infty}),
set breakpoints = c(0,1,12), or breakpoints = c(1,12).}

\item{init_vals}{an optional vector of initial values at which to start the optimizer for the
parameters. Must specify the appropriate number of parameters for the given distribution /
number of periods.}

\item{etsp_c}{if "etsp" distribution is used, the value to set the c parameter equal to. Defaults to 0.}

\item{force_nonincreasing}{a boolean, specifying (for the log-logistic and piecewise exponential
distributions) whether to force the hazard to be non-increasing across age. Defaults to FALSE. This option
is not available for other distributions.}
}
\value{
A list containing:
\itemize{
\item result: a dataframe of summarized results
order arrange(time). Each matrix will have rows arranged in order arrange(region).
\item optim: the output from \code{optim}
\item grad: the gradient evaluated at the MLE
\item variance: the finite population variance-covariance matrix
\item design: the survey design object
\item runtime: runtime for likelihood optimization
\item initial_values: initial values for parameters used in likelihood maximization
}
}
\description{
Fits a parametric survival model to DHS data (formatted from \code{format_dhs}) to get survey-weighted
parameter estimates and finite population variances for parametric survival curves in each time period.
}
\details{

}
\references{
Stacy, E. W. (1962). A generalization of the gamma
distribution.  Annals of Mathematical Statistics 33:1187-92.

Prentice, R. L. (1974). A log gamma model and its maximum likelihood
estimation. Biometrika 61(3):539-544.

Scholey, J (2019). The Age-Trajectory of Infant Mortality in the United States:
Parametric Models and Generative Mechanisms. Annual meeting of the Population
Association of America, Austin, TX. 2019.
}
\author{
Taylor Okonek
}
