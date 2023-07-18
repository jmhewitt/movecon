#include "ProjectedLocationLikelihood.h"
#include "DomainSearch.h"

/**
 * Create a family of location observation distributions from ellipse vectors
*/
std::vector<ProjectedLocationLikelihood> LocationDistributionFamily(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations
) {
    std::vector<ProjectedLocationLikelihood> family;
    family.reserve(eastings.size());
    
    auto eastings_it = eastings.begin();
    auto northings_it = northings.begin();
    auto semi_majors_it = semi_majors.begin();
    auto semi_minors_it = semi_minors.begin();
    auto orientations_it = orientations.begin();

    auto eastings_end = eastings.end();
    for(; eastings_it != eastings_end; ++eastings_it) {
        family.emplace_back(
            ProjectedLocationLikelihood::from_ellipse(
                *eastings_it, *(northings_it++), *(semi_majors_it++), 
                *(semi_minors_it++), *(orientations_it++)
            )
        );
    }

    return family;
}

/**
 * Create a family of location observation distributions from GPS vectors
*/
std::vector<ProjectedLocationLikelihood> LocationDistributionFamilyFromGPS(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> hdops, double uere
) {
    std::vector<ProjectedLocationLikelihood> family;
    family.reserve(eastings.size());
    
    auto eastings_it = eastings.begin();
    auto northings_it = northings.begin();
    auto hdops_it = hdops.begin();

    auto eastings_end = eastings.end();
    for(; eastings_it != eastings_end; ++eastings_it) {
        family.emplace_back(
            ProjectedLocationLikelihood::from_hdop_uere(
                *eastings_it, *(northings_it++), *(hdops_it++), uere
            )
        );
    }

    return family;
}

/**
 * Sample states from a Gaussian distribution constrained to a spatial domain,
 * and with last movement directions uniformly sampled
 * 
*/
// [[Rcpp::export]]
Rcpp::List sample_gaussian_states(
    Rcpp::XPtr<RookDirectionalStatespaceSearch> statespace_search, 
    double easting,
    double northing,
    double semi_major,
    double semi_minor,
    double orientation,
    std::size_t n
) {

    typedef RookDirectionalStatespace::StateType StateType;

    // initialize sampler
    ProjectedLocationLikelihood sampler = 
        ProjectedLocationLikelihood::from_ellipse(
            easting, northing, semi_major, semi_minor, orientation
        );

    //
    // initialize output 
    //

    std::vector<StateType*> * states = new std::vector<StateType*>();
    states->reserve(n);

    Rcpp::List states_formatted;

    //
    // sample states
    //

    for(std::size_t i = 0; i < n; ++i) {
        // draw gaussian coordinates
        double r_easting, r_northing;
        sampler.sample(r_easting, r_northing);
        // map to grid
        Location * r_location = statespace_search->map_location(
            r_easting, r_northing
        );
        // randomly select a state at location
        std::set<StateType*> r_states_set = 
            statespace_search->states_by_location[r_location];
        std::size_t r_ind = std::floor(R::runif(0, r_states_set.size()));
        auto set_iter = r_states_set.begin();
        for(std::size_t j = 0; j < r_ind; ++j) 
            ++set_iter;
        StateType* r_state = *set_iter;
        // export sample
        states->push_back(r_state);
        // states_formatted.push_back(format_state(*r_state));
    }

    // package results
    Rcpp::XPtr<std::vector<StateType*>> state_ptr(states, true);
    return Rcpp::List::create(
        Rcpp::Named("states") = states_formatted,
        Rcpp::Named("states_cpp") = state_ptr
    );
}
