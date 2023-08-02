#include "ParticleGillespie.h"
#include "Domain.h"
#include "Tx.h"
#include "Directions.h"

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

/**
 * Forward-simulate movement on a statespace using Gillespie
*/
// [[Rcpp::export]]
Rcpp::List Test__Particle_Gillespie_Steps(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    std::string last_movement_direction,
    std::size_t easting_ind, 
    std::size_t northing_ind,
    double directional_persistence,
    Eigen::VectorXd beta,
    std::vector<double> times
) {

    // get starting state
    typedef RookDirectionalStatespace::StateKey StateKey;
    typedef RookDirectionalStatespace::StateType StateType;
    StateType & state = statespace->states.at(
        StateKey(
            stringToDirection(last_movement_direction), 
            easting_ind, 
            northing_ind
        )
    );

    // construct transition rate evaluator
    typedef location_based_movement<StateType, Eigen::VectorXd> 
        particle_transition_rate;
    particle_transition_rate location_based_rate(beta);

    // construct transition probability evaluator
    typedef directional_transition_probabilities<
        StateType, CardinalDirectionOrientations
    > particle_transition_probability;
    particle_transition_probability transition_prob(directional_persistence);

    // build a particle at the state
    ParticleGillespie<
        StateType, 
        particle_transition_rate,
        particle_transition_probability
    > particle(location_based_rate, transition_prob);
    particle.state = &state;

    // initialize output
    Rcpp::List path;
    path.push_back(format_state(state));

    // run forward simulation
    for(auto it = (times.begin()+1); it != times.end(); ++it) {
        particle.step(*(it-1), *it);
        path.push_back(format_state(*particle.state));
    }

    // return simulated path
    return path;
}
