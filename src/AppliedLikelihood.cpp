#include "AppliedLikelihood.h"

/**
 * Create a family of location observation distributions from input vectors in 
 * which location observations are not necessarily available at all discrete 
 * timepoints.
 * 
 * Flat likelihoods will be inserted into the family at times without 
 * location observations.
 * 
 * The return type is an std::vector of std::unique_ptr likelihood objects to 
 * mitigate object slicing that would occur by combining different types of 
 * objects that have the same interface for member functions but different 
 * implementations (i.e., to model different likelihood contributions).
 * 
 * @param t vector of discrete time indices (starting at 0) at which 
 *   observations are available
 * @param nt total number of discrete timepoints
 * 
*/
std::vector<std::unique_ptr<AppliedLikelihood>> AppliedLikelihoodFamily(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations, std::vector<std::size_t> t,
    std::size_t nt
) {
    std::vector<std::unique_ptr<AppliedLikelihood>> family;
    family.reserve(eastings.size());
    
    auto t_it = t.begin();
    auto eastings_it = eastings.begin();
    auto northings_it = northings.begin();
    auto semi_majors_it = semi_majors.begin();
    auto semi_minors_it = semi_minors.begin();
    auto orientations_it = orientations.begin();

    for(std::size_t ind = 0; ind < nt; ++ind) {
        if(ind == *t_it) {
            ++t_it;
            family.emplace_back(
                new AppliedLocationLikelihood(
                    AppliedLocationLikelihood::from_ellipse(
                        *(eastings_it++), *(northings_it++), 
                        *(semi_majors_it++), *(semi_minors_it++), 
                        *(orientations_it++)
                    )
                )
            );
        } else {
            family.emplace_back(new AppliedFlatLikelihood());
        }
    }

    return family;
}

std::vector<std::unique_ptr<AppliedLikelihood>> AppliedLikelihoodFamilyFromGPS(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> hdops, double uere, std::vector<std::size_t> t,
    std::size_t nt
) {
    std::vector<std::unique_ptr<AppliedLikelihood>> family;
    family.reserve(eastings.size());
    
    auto t_it = t.begin();
    auto eastings_it = eastings.begin();
    auto northings_it = northings.begin();
    auto hdops_it = hdops.begin();

    for(std::size_t ind = 0; ind < nt; ++ind) {
        if(ind == *t_it) {
            ++t_it;
            family.emplace_back(
                new AppliedLocationLikelihood(
                    AppliedLocationLikelihood::from_hdop_uere(
                        *(eastings_it++), *(northings_it++), *(hdops_it++), uere
                    )
                )
            );
        } else {
            family.emplace_back(new AppliedFlatLikelihood());
        }
    }

    return family;
}

/**
 * Demonstrate that we can create and call likelihoods of mixed types
*/
// [[Rcpp::export]]
double Test__AppliedLikelihoodFamily(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations, std::vector<std::size_t> t,
    std::size_t nt,
    Rcpp::XPtr<
        std::vector<RookDirectionalStatespace::StateType*>
    > states
) {

    std::vector<std::unique_ptr<AppliedLikelihood>> f = AppliedLikelihoodFamily(
        eastings, northings, semi_majors, semi_minors, orientations, t, nt
    );
    
    return f[0]->dstate(**states->begin()) + f[1]->dstate(**states->begin());
}
