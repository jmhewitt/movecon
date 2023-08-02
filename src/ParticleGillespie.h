#ifndef MOVECON_PARTICLE_GILLESPIE_H
#define MOVECON_PARTICLE_GILLESPIE_H

#include <Rcpp.h>

template<
    typename StateType, 
    // Type that can evaluate Hewitt et. al. (2023) eq. 14
    typename transition_rate_evaluator,
    // Type that can evaluate Hewitt et. al. (2023) eq. 15
    typename transition_probability_evaluator
>
struct ParticleGillespie {

    private:

        transition_rate_evaluator* m_rate_evaluator;
        transition_probability_evaluator* m_probability_evaluator;

    public:

        StateType* state; 

        ParticleGillespie(
            transition_rate_evaluator & rate_evaluator,
            transition_probability_evaluator & probability_evaluator
        ) : m_rate_evaluator(&rate_evaluator), 
            m_probability_evaluator(&probability_evaluator) { }

        /**
         * Forward simulation via Gillespie algorithm.
        */
        void step(double t, double tnext) {
            // initial time increment
            t += R::rexp(1 / m_rate_evaluator->transition_rate(*state));
            // transition to neighbors while able (i.e., before tnext)
            while(t < tnext) {
                // get transition probabilities
                const double * mass = 
                    m_probability_evaluator->probabilities(*state).data();
                // transition to random neighbor
                double p = R::runif(0, 1);
                double cumulative_mass = 0;
                for(auto destination : state->to) {
                    // aggregate transition mass from neighbor
                    cumulative_mass += *(mass++);
                    if(cumulative_mass > p) {
                        // transition to neighbor
                        state = destination;
                        break;
                    }
                }
                // increment time
                t += R::rexp(1 / m_rate_evaluator->transition_rate(*state));
            }
        } // step function

}; 

#endif
