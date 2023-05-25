#ifndef MOVECON_PARTICLE_H
#define MOVECON_PARTICLE_H

#include <Rcpp.h>

template<
    typename StateType, 
    // Type that can evaluate Hewitt et. al. (2023) eq. 14
    typename transition_rate_evaluator,
    // Type that can evaluate Hewitt et. al. (2023) eq. 15
    typename transition_probability_evaluator
>
struct Particle {

    private:

        transition_rate_evaluator* m_rate_evaluator;
        transition_probability_evaluator* m_probability_evaluator;

    public:

        StateType* state; 

        Particle(
            transition_rate_evaluator & rate_evaluator,
            transition_probability_evaluator & probability_evaluator
        ) : m_rate_evaluator(&rate_evaluator), 
            m_probability_evaluator(&probability_evaluator) { }

        /**
         * forward-simulation using discretized transition distribution
        */
        void step() {

            double uniformized_rate = m_rate_evaluator->transition_rate(*state);
            
            // test for self-transition, then move
            if(R::runif(0, 1) < 1 - uniformized_rate) {
                // self-transition, do nothing
            } else {
                
                // get transition probabilities
                auto probs = m_probability_evaluator->probabilities(*state);
                
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
