#include "Directions.h"

/**
 * Convert a direction string to a CTDS domain object type
*/
CardinalDirection stringToDirection(const std::string & direction) {
    if(direction.compare("north") == 0) {
        return CardinalDirection::north;
    } else if(direction.compare("east") == 0) {
        return CardinalDirection::east;
    } else if(direction.compare("south") == 0) {
        return CardinalDirection::south;
    } else if(direction.compare("west") == 0) {
        return CardinalDirection::west;
    }
    Rcpp::stop("Argument direction contains an invalid value");
}

/**
 * Format a CTDS domain object direction type for printing
*/
std::string directionToString(const CardinalDirection & direction) {
    switch(direction) {
        case CardinalDirection::north:
            return "north";
        case CardinalDirection::east:
            return "east";
        case CardinalDirection::south:
            return "south";
        case CardinalDirection::west:
            return "west";
    }
    Rcpp::stop("String conversion is not defined for argument direction");
}

// [[Rcpp::export]]
double TestDirectionalCovariate(std::string x, std::string y) {
    return CardinalDirectionOrientations().directional_persistence_covariate(
        stringToDirection(x), stringToDirection(y)
    );
}