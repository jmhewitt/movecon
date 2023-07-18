#include "ParticleFilter.h"

#include "Particle.h"
#include "Domain.h"
#include "Tx.h"
#include "Directions.h"
#include "AppliedLikelihood.h"

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

Rcpp::List run_particle_filter(
    /* likelihood components */
    std::vector<std::unique_ptr<AppliedLikelihood>> & likelihood_seq,
    /* filter components */
    Rcpp::XPtr<RookDirectionalStatespace> & statespace,
    Rcpp::XPtr<
        std::vector<RookDirectionalStatespace::StateType*>
    > & initial_latent_state_sample,
    /* model parameters */
    double directional_persistence, Eigen::VectorXd & beta, double delta
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

    typedef RookDirectionalStatespace StatespaceType;

    typedef StatespaceType::StateKey StateKey;
    typedef StatespaceType::StateType StateType;

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

    typedef std::vector<std::unique_ptr<AppliedLikelihood>> LikelihoodSeqType;

    //
    // build particles
    //

    // construct transition rate evaluator
    base_transition_rate location_based_rate(beta);
    uniformized_transition_rate uniformized_rate(
        &location_based_rate, delta
    );
    particle_transition_rate transition_rate(uniformized_rate);

    // construct transition probability evaluator
    directional_probabilities directional_probs(directional_persistence);
    particle_transition_probability transition_prob(directional_probs);

    // initialize particle container
    std::vector<ParticleType> particles;
    particles.reserve(initial_latent_state_sample->size());

    // build particles for the initial states
    ParticleType particle(transition_rate, transition_prob);
    auto state = initial_latent_state_sample->begin();
    auto state_end = initial_latent_state_sample->end();
    for(; state != state_end; ++state) {
        particle.state = *state;
        particles.push_back(particle);
    }

    //
    // build proposal distributions
    //

    ProposalSeqType proposal_seq = ConstantStepFamily<ParticleType>(
        likelihood_seq.size(), 1
    );

    //
    // build particle filter
    //

    BootstrapParticleFilter<
        ParticleType, 
        ProposalSeqType, 
        LikelihoodSeqType,
        FilterObserver<ParticleType>
    > 
    pf(particles);

    pf.proposal_distributions = &proposal_seq;
    pf.likelihoods = &likelihood_seq;

    // raw storage for filtering distributions
    FilterObserver<ParticleType> filtering_distributions;

    // run particle filter
    double ll = pf.marginal_ll(filtering_distributions);

    //
    // export filtering distributions as an array
    //

    Rcpp::NumericVector filtering_locations(
        Rcpp::Dimension(
            2, // coordinates
            particles.size(), // particles
            filtering_distributions.particle_distributions.size() // dist'ns.
        )
    );

    // export filtering distributions as coordinates
    double * filtering_loc = filtering_locations.begin();
    auto distribution = filtering_distributions.particle_distributions.begin();
    auto dist_end = filtering_distributions.particle_distributions.end();
    for(; distribution != dist_end; ++distribution) {
        // loop over particles within each distribution
        auto particle = distribution->begin();
        auto particle_end = distribution->end();
        for(; particle != particle_end; ++particle) {
            // transfer coordinates
            *(filtering_loc++) = particle->state->properties.location->easting;
            *(filtering_loc++) = particle->state->properties.location->northing;
        } // particle
    } // distribution


    // package results
    return Rcpp::List::create(
        Rcpp::Named("ll") = ll,
        Rcpp::Named("filtering_distributions") = filtering_locations
    );
}

// [[Rcpp::export]]
Rcpp::List Test__Particle_Filter_Likelihood(
    /* likelihood components */
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations,
    std::vector<std::size_t> t,
    std::size_t nt,
    /* filter components */
    Rcpp::XPtr<RookDirectionalStatespace> statespace,
    Rcpp::XPtr<
        std::vector<RookDirectionalStatespace::StateType*>
    > initial_latent_state_sample,
    /* model parameters */
    double directional_persistence, Eigen::VectorXd beta, double delta
) {
    
    typedef std::vector<std::unique_ptr<AppliedLikelihood>> LikelihoodSeqType;

    LikelihoodSeqType likelihood_seq = AppliedLikelihoodFamily(
        eastings, northings, semi_majors, semi_minors, orientations, t, nt
    );

    return run_particle_filter(
        likelihood_seq, statespace, initial_latent_state_sample,
        directional_persistence, beta, delta
    );
}
