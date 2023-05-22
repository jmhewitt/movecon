/**
 * Ensure user-defined C++ types are accessible within RcppExports.cpp when 
 * custom types are used as arguments for Rcpp::export functions (i.e., for 
 * working with Rcpp::XPtr objects)
 * 
 * source: Section 2.5 in Rcpp attributes description at 
 *   https://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-attributes.pdf
*/ 
#include "Domain.h"
