# pssst

Functions to fit parametric survival models to survey data to estimate finite-population, synthetic child mortality across multiple time periods.

## Overview

**pssst** is an open-source R package containing functions to fit parametric survival models to survey data to estimate finite-population, synthetic child mortality across multiple time periods. 

The package currently supports the following types of survival data:

* Interval-censored
* Right-censored
* Exactly observed

The package currently supports the following parametric distributions:

* Exponential
* Weibull
* Piecewise Exponential

## Installation

The current development version can be installed using `devtools::install_github()`:

```R
devtools::install_github(repo="taylorokonek/pssst")
```


