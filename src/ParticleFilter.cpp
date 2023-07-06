#include "ParticleFilter.h"

#include "Particle.h"
#include "Domain.h"
#include "Tx.h"
#include "Directions.h"
#include "log_add.h"
#include "ProjectedLocationLikelihood.h"

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
double Test__Particle_Filter_Likelihood (

    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations,

    Rcpp::XPtr<RookDirectionalStatespace> statespace,
    std::string last_movement_direction,
    std::size_t easting_ind, 
    std::size_t northing_ind,
    std::size_t nparticles,

    double directional_persistence, Eigen::VectorXd beta, double delta

) {

    // reset cached state values
    auto end = statespace->states.end();
    for(auto state = statespace->states.begin(); state != end; ++state) {
        state->second.to_rate = -1;
        state->second.to_probabilities.resize(0);
    }

    //
    // configurations
    //

    typedef RookDirectionalStatespace::StateKey StateKey;
    typedef RookDirectionalStatespace::StateType StateType;

    typedef location_based_movement<StateType, Eigen::VectorXd> 
        base_transition_rate;
    typedef uniformized_rate_evaluator<StateType, base_transition_rate> 
        uniformized_transition_rate;
    typedef state_cache_rate_evaluator<StateType,  uniformized_transition_rate> 
        particle_transition_rate;

    typedef directional_transition_probabilities<
        StateType, CardinalDirectionOrientations
    > directional_probabilities;
    typedef state_cache_transition_probability_evaluator<
        StateType, directional_probabilities
    > particle_transition_probability;

    typedef Particle<
        StateType, 
        particle_transition_rate,
        particle_transition_probability
    > ParticleType;

    typedef std::vector<NStepProposal<ParticleType>> ProposalSeqType;

    typedef std::vector<
        ProjectedLocationLikelihood<ParticleType>
    > LikelihoodSeqType;

    //
    // build particles
    //

    // get starting state
    StateType & state = statespace->states.at(
        StateKey(
            stringToDirection(last_movement_direction), 
            easting_ind, 
            northing_ind
        )
    );

    // construct transition rate evaluator
    base_transition_rate location_based_rate(beta);
    uniformized_transition_rate uniformized_rate(
        &location_based_rate, delta
    );
    particle_transition_rate transition_rate(uniformized_rate);

    // construct transition probability evaluator
    directional_probabilities directional_probs(directional_persistence);
    particle_transition_probability transition_prob(directional_probs);

    // build a particle at the state
    ParticleType particle(transition_rate, transition_prob);
    particle.state = &state;

    // TODO: better initialization for particles
    std::vector<ParticleType> particles(nparticles, particle);

    //
    // build likelihoods
    //

    LikelihoodSeqType likelihood_seq = LocationDistributionFamily<ParticleType>(
        eastings, northings, semi_majors, semi_minors, orientations
    );

    //
    // build proposal distributions
    //

    ProposalSeqType proposal_seq = ConstantStepFamily<ParticleType>(
        eastings.size(), 1
    );

    //
    // build and run particle filter
    //

    BootstrapParticleFilter<ParticleType, ProposalSeqType, LikelihoodSeqType> 
    pf(particles);

    pf.proposal_distributions = &proposal_seq;
    pf.likelihoods = &likelihood_seq;

    return pf.marginal_ll();
}