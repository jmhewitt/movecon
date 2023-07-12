#include "AppliedLikelihood.h"

/**
 * Demonstrate that we can create and call likelihoods of mixed types
*/
// [[Rcpp::export]]
double Test__AppliedLikelihoodFamily(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations, std::vector<std::size_t> t,
    Rcpp::XPtr<
        std::vector<RookDirectionalStatespace::StateType*>
    > states
) {

    std::vector<std::unique_ptr<AppliedLikelihood>> f = AppliedLikelihoodFamily(
        eastings, northings, semi_majors, semi_minors, orientations
    );
    
    return f[0]->dstate(**states->begin()) + f[1]->dstate(**states->begin());
}
