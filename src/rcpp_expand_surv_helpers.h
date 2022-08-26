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

