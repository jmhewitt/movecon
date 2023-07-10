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

// [[Rcpp::export]]
Rcpp::List Test__States_At_Location(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    double easting,
    double northing
) {

    typedef RookDirectionalStatespace::StateType StateType;

    StatespaceSearch<RookDirectionalStatespace> search(*statespace);

    Location * location = search.map_location(easting, northing);
    std::set<StateType*> states = search.states_by_location[location];

    Rcpp::List res;

    auto state = states.begin();
    auto state_end = states.end();
    for(; state != state_end; ++state) {
        res.push_back(format_state(**state));
    }

    return res;
}
