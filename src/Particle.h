#ifndef MOVECON_PARTICLE_H
#define MOVECON_PARTICLE_H

#include <Rcpp.h>

template<
    typename StateType, 
    // Type that can evaluate Hewitt et. al. (2023) eq. 14
    typename transition_rate_builder,
    // Type that can evaluate Hewitt et. al. (2023) eq. 15
    typename transition_probability_builder
>
struct Particle {

    StateType* state; 

    /**
     * forward-simulation using discretized transition distribution
    */
    template<typename VectorType>
    void step(
        const VectorType & beta, 
        const double delta, 
        const double directional_persistence
    ) {

        double uniformized_rate = delta * 
            transition_rate_builder::transition_rate(*state, beta);
        
        // test for self-transition, then move
        if(R::runif(0, 1) < 1 - uniformized_rate) {
            // self-transition, do nothing
        } else {
            
            // get transition probabilities
            auto probs = transition_probability_builder::probabilities(
                *state, directional_persistence
            );
            
            // transition to random neighbor
            double p = R::runif(0, 1);
            double cumulative_mass = 0;
            double * mass = probs.data();
            for(auto destination : state->to) {
                // aggregate transition mass from neighbor
                cumulative_mass += *(mass++);
                if(cumulative_mass > p) {
                    // transition to neighbor
                    state = destination;
                    break;
                }
            }
        } // transition logic
    } // step function

};

#endif
