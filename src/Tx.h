/**
 * Objects that define a local transition distribution
*/

#ifndef MOVECON_TX_H
#define MOVECON_TX_H

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

template<typename State, typename DirectionalPersistence>
struct directional_transition_probabilities {

    /**
     * Evaluate Hewitt et. al. (2023) eq. 15, specialized for directional 
     * persistence as the only directional driver of movement.
     * 
     * Returns the probability that each of the adjacent states will be visited 
     * when transitioning away from the state argument
     * 
     * @param state State object that has information about directional 
     *   persistence and its neighbors
     * @param directional_persistence scalar that indicates the strength of 
     *   directional persistence for a transition.  Use 0 for random walks, 
     *   which lack directional persistence.
     * 
    */
    static Eigen::VectorXd probabilities(
        const State & state, const double directional_persistence
    ) {

        // initialize output
        Eigen::VectorXd p(state.to.size());
        double* pIt = p.data();

        // aggregate probability transition mass for states that can be reached
        auto stateEnd = state.to.end();
        for(auto stateIt = state.to.begin(); stateIt != stateEnd; ++stateIt) {
            *(pIt++) = std::exp(
                directional_persistence * 
                DirectionalPersistence::directional_persistence_covariate(
                    state.properties.last_movement_direction,
                    (*stateIt)->properties.last_movement_direction
                )
            );
        }

        // standardize transition distribution
        return p / p.sum();
    }

};

#endif
