/**
 * Objects that define a local transition distribution
*/

#ifndef MOVECON_TX_H
#define MOVECON_TX_H

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

template<typename State, typename DirectionalPersistence>
class directional_transition_probabilities {

    private:

        const double directional_persistence;

    public:

        /**
         * @param persistence scalar that indicates the strength of 
         *   directional persistence for a transition.  Use 0 for random walks, 
         *   which lack directional persistence.
        */
        directional_transition_probabilities(double persistence) : 
            directional_persistence(persistence) { }

        /**
         * Evaluate Hewitt et. al. (2023) eq. 15, specialized for directional 
         * persistence as the only directional driver of movement.
         * 
         * Returns the probability that each of the adjacent states will be  
         * visited when transitioning away from the state argument
         * 
         * @param state State object that has information about directional 
         *   persistence and its neighbors
         * 
        */
        Eigen::VectorXd probabilities(const State & state) {

            // initialize output
            Eigen::VectorXd p(state.to.size());
            double* pIt = p.data();

            // aggregate probability tx. mass for states that can be reached
            auto stateIt = state.to.begin();
            auto stateEnd = state.to.end();
            for(; stateIt != stateEnd; ++stateIt) {
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

template<typename State, typename VectorType>
struct location_based_movement {

    const VectorType & beta;

    location_based_movement(const VectorType & x) : beta(x) { }

    /**
     * Evaluate Hewitt et. al. (2023) eq. 14, which specifies a total transition
     * rate using location-based covariates
    */
   double transition_rate(const State & state) {
        return std::exp(beta.dot(state.properties.location->x));
   }

};

/**
 * Scale transition rate from a transition_rate_evaluator object by a constant
*/
template<typename State, typename transition_rate_evaluator>
class uniformized_rate_evaluator {

    private:

        transition_rate_evaluator* m_evaluator;
        double m_scale;

    public: 

        uniformized_rate_evaluator(
            transition_rate_evaluator* evaluator, double scale
        ) : m_evaluator(evaluator), m_scale(scale) { }

        double transition_rate(const State & state) {
            return m_scale * m_evaluator->transition_rate(state);
        }
};

#endif
