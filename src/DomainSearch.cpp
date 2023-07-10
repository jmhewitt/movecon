#include "DomainSearch.h"

// [[Rcpp::export]]
Rcpp::List Test__Map_Location(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    double easting,
    double northing
) {
    StatespaceSearch<RookDirectionalStatespace> search(*statespace);
    return format_location(*search.map_location(easting, northing));
}