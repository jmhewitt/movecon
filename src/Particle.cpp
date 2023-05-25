#include "Particle.h"
#include "Domain.h"
#include "Tx.h"
#include "Directions.h"

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

/**
 * Forward-simulate movement on a statespace
*/
// [[Rcpp::export]]
Rcpp::List Test__Particle_Steps(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    std::string last_movement_direction,
    std::size_t easting_ind, 
    std::size_t northing_ind,
    double directional_persistence,
    Eigen::VectorXd beta,
    double delta,
    std::size_t nsteps
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

    // build a particle at the state
    Particle<
        StateType, 
        location_based_movement<StateType>,
        directional_transition_probabilities< 
            StateType, CardinalDirectionOrientations
        >
    > particle;
    particle.state = &state;

    // initialize output
    Rcpp::List path;
    path.push_back(format_state(state));

    // run forward simulation
    for(std::size_t i = 0; i < nsteps; ++i) {
        particle.step(beta, delta, directional_persistence);
        path.push_back(format_state(*particle.state));
    }

    // return simulated path
    return path;
}
