#include "Tx.h"
#include "Domain.h"
#include "CacheDecorator.h"

/**
 * Compute transition probabilities for a given state in a statespace
*/
// [[Rcpp::export]]
Eigen::VectorXd Test__Directional_Transition_Probabilities(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    std::string last_movement_direction,
    std::size_t easting_ind, 
    std::size_t northing_ind,
    double directional_persistence
) {
    
    typedef RookDirectionalStatespace::StateKey StateKey;
    typedef RookDirectionalStatespace::StateType StateType;
    StateType & state = statespace->states.at(
        StateKey(
            stringToDirection(last_movement_direction), 
            easting_ind, 
            northing_ind
        )
    );

    directional_transition_probabilities<
        StateType, CardinalDirectionOrientations
    > probability_builder(directional_persistence);

    return probability_builder.probabilities(state);
}

// /**
//  * Repeatedly compute transition probabilities for a given state in a 
//  * statespace, using caching to avoid redundant computations
// */
// // [[Rcpp::export]]
// Eigen::VectorXd Test__Directional_Transition_Probability_Cache(
//     Rcpp::XPtr<RookDirectionalStatespace> statespace, 
//     std::string last_movement_direction,
//     std::size_t easting_ind, 
//     std::size_t northing_ind,
//     double directional_persistence,
//     std::size_t reps
// ) {
    
//     typedef RookDirectionalStatespace::StateKey StateKey;
//     typedef RookDirectionalStatespace::StateType StateType;
//     StateType & state = statespace->states.at(
//         StateKey(
//             stringToDirection(last_movement_direction), 
//             easting_ind, 
//             northing_ind
//         )
//     );

//     typedef directional_transition_probabilities<
//         StateType, CardinalDirectionOrientations
//     > probability_builder;

//     probability_builder pb(directional_persistence);

//     CacheDecorator<Eigen::VectorXd, StateType, double> cd(
//         &pb::probabilities
//     );

//     for(std::size_t i = 0; i < reps; ++i) {
//         cd(state);
//     }

//     return cd(state);
// }

// /**
//  * Repeatedly compute transition probabilities for a given state in a 
//  * statespace, meant to facilitate comparison with performance using caching
// */
// // [[Rcpp::export]]
// Eigen::VectorXd Test__Directional_Transition_Probability_No_Cache(
//     Rcpp::XPtr<RookDirectionalStatespace> statespace, 
//     std::string last_movement_direction,
//     std::size_t easting_ind, 
//     std::size_t northing_ind,
//     double directional_persistence,
//     std::size_t reps
// ) {
    
//     typedef RookDirectionalStatespace::StateKey StateKey;
//     typedef RookDirectionalStatespace::StateType StateType;
//     StateType & state = statespace->states.at(
//         StateKey(
//             stringToDirection(last_movement_direction), 
//             easting_ind, 
//             northing_ind
//         )
//     );

//     typedef directional_transition_probabilities<
//         StateType, CardinalDirectionOrientations
//     > probability_builder;

//     probability_builder pb(directional_persistence);


//     CacheDecorator<Eigen::VectorXd, StateType, double> cd(
//         std::mem_fn(&probability_builder::probabilities)
//     );

//     for(std::size_t i = 0; i < reps; ++i) {
//         pb.probabilities(state);
//     }

//     return pb.probabilities(state);
// }

/**
 * Compute transition rate for a given state in a statespace
*/
// [[Rcpp::export]]
double Test__Location_Based_Movement_Transition_Rate(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    std::string last_movement_direction,
    std::size_t easting_ind, 
    std::size_t northing_ind,
    Eigen::VectorXd beta
) {

    typedef RookDirectionalStatespace::StateKey StateKey;
    typedef RookDirectionalStatespace::StateType StateType;
    StateType & state = statespace->states.at(
        StateKey(
            stringToDirection(last_movement_direction), 
            easting_ind, 
            northing_ind
        )
    );

    location_based_movement<
        StateType, Eigen::VectorXd
    > transition_rate_builder(beta);

    return transition_rate_builder.transition_rate(state);
}
