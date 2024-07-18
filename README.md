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
* Generalized Gamma
* Lognormal
* Gompertz
* Log-logistic

## Installation

The current development version can be installed using `devtools::install_github()`:

```R
devtools::install_github(repo="taylorokonek/pssst")
```

### Common Installation Error

If you receive an error that contains the following message...

```R
ld: warning: directory not found for option '-L/opt/gfortran/lib/gcc/aarch64-apple-darwin20.0/12.2.0'
ld: warning: directory not found for option '-L/opt/gfortran/lib'
ld: library not found for -lgfortran
```

... you likely have a Mac with an M1 chip, and it is your first time installing a package that depends on c++ code! The fix for this is *typically* the following (steps taking from <a href="https://github.com/Rdatatable/data.table/wiki/Installation" >this</a> post):


1. Create a `~/.R/Makevars` file on your machine, if one does not already exist

2. Open the Makevars file in your favorite text editor, and add the following lines to it:

```R
LOC = /usr/local/gfortran
CC=$(LOC)/bin/gcc -fopenmp
CXX=$(LOC)/bin/g++ -fopenmp
CXX11 = $(LOC)/bin/g++ -fopenmp # for fst package
CFLAGS=-g -O3 -Wall -pedantic -std=gnu99 -mtune=native -pipe
CXXFLAGS=-g -O3 -Wall -pedantic -std=c++11 -mtune=native -pipe
LDFLAGS=-L$(LOC)/lib -Wl,-rpath,$(LOC)/lib
CPPFLAGS=-I$(LOC)/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include
```

3. Save the file, and try re-installing **pssst**. If it still doesn't work, please create an Issue on this github page so we can address this!






