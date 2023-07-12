#ifndef MOVECON_APPLIED_LIKELIHOOD_H
#define MOVECON_APPLIED_LIKELIHOOD_H

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

#include "Tx.h"
#include "Particle.h"
#include "Directions.h"
#include "Domain.h"

#include "ProjectedLocationLikelihood.h"

struct AppliedLikelihood {

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

    virtual double dparticle(const ParticleType & particle) = 0;

    virtual double dstate(const StateType & state) = 0;

};

struct AppliedFlatLikelihood : public AppliedLikelihood {
    double dparticle(const ParticleType & particle) {
        return 0;
    }
    double dstate(const StateType & state) {
        return 0;
    }
};

struct AppliedLocationLikelihood : public AppliedLikelihood {

    private:

        ProjectedLocationLikelihood likelihood_impl;

    public:

        AppliedLocationLikelihood(
            double easting, double northing, double semi_major, 
            double semi_minor, double orientation
        ) : likelihood_impl(
                easting, northing, semi_major, semi_minor, orientation
            ) { }

        double dparticle(const ParticleType & particle) {
            return likelihood_impl.dparticle(particle);
        }

        double dstate(const StateType & state) {
            return likelihood_impl.dstate(state);
        }

};

std::vector<std::unique_ptr<AppliedLikelihood>> AppliedLikelihoodFamily(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations, std::vector<std::size_t> t,
    std::size_t nt
);

#endif