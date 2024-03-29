#include "DomainSearch.h"

// [[Rcpp::export]]
Rcpp::XPtr<RookDirectionalStatespaceSearch> build_statespace_search(
    Rcpp::XPtr<RookDirectionalStatespace> statespace
) {
    RookDirectionalStatespaceSearch * search = 
        new RookDirectionalStatespaceSearch(*statespace);
    return Rcpp::XPtr<RookDirectionalStatespaceSearch>(search, true);
}

// [[Rcpp::export]]
Rcpp::List nearest_location_in_domain(
    Rcpp::XPtr<RookDirectionalStatespaceSearch> statespace_search, 
    double easting,
    double northing
) {
    return format_location(*statespace_search->map_location(easting, northing));
}

// [[Rcpp::export]]
Rcpp::List states_at_nearest_location_in_domain(
    Rcpp::XPtr<RookDirectionalStatespaceSearch> statespace_search, 
    double easting,
    double northing
) {

    typedef RookDirectionalStatespace::StateType StateType;

    Location * location = statespace_search->map_location(easting, northing);
    std::set<StateType*> states = 
        statespace_search->states_by_location[location];

    Rcpp::List res;

    auto state = states.begin();
    auto state_end = states.end();
    for(; state != state_end; ++state) {
        res.push_back(format_state(**state));
    }

    return res;
}
