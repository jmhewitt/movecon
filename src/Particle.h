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
            } // transition logic
        } // step function

        /**
         * n-steps of forward-simulation using discretized tx. distribution
        */
        void step(std::size_t n) {
            for(std::size_t i = 0; i < n; ++i)
                step();
        }

};

/**
 * Wrapper to call a particle's member step function nstep times.
*/
template<typename Particle>
class NStepProposal {
    
    private: 

        std::size_t nsteps;

    public:

        NStepProposal(std::size_t n) : nsteps(n) { }

        void propose(Particle & particle) {
            particle.step(nsteps);
        }

};

/**
 * Create a family of identical proposal distributions
 * 
 * Although this is somewhat inefficient with memory, it is simple to implement
 * 
 * @param size Number of proposal distribution copies to put in the family
 * @param nsteps Number of simulation steps for each proposal distribution
*/
template<typename Particle>
std::vector<NStepProposal<Particle>> ConstantStepFamily(
    std::size_t size, std::size_t nsteps
) {
    NStepProposal<Particle> proposal(nsteps);
    std::vector<NStepProposal<Particle>> family(size, proposal);
    return family;
}

/**
 * Create a family of proposal distributions from a vector of step sizes
 * 
 * @param steps Vector specifying number of steps for each proposal distribution
*/
template<typename Particle, typename stepType>
std::vector<NStepProposal<Particle>> DiscretizedTimestepFamily(
    std::vector<std::size_t> steps
) {
    typedef NStepProposal<Particle> ProposalDistribution;
    std::vector<ProposalDistribution> family(steps.size());
    auto end = steps.end();
    for(auto it = steps.begin(); it != end; ++it) {
        ProposalDistribution tmp(*it);
        family.emplace_back(tmp);
    }
    return family;
}

#endif
