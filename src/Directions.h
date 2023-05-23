#ifndef MOVECON_DIRECTIONS_H
#define MOVECON_DIRECTIONS_H

#include <Rcpp.h>

enum CardinalDirection { north = 0, east = 1, south = 2, west = 3 };

std::string directionToString(const CardinalDirection & direction);
CardinalDirection stringToDirection(const std::string & direction);

class CardinalDirectionOrientations {

    private:

        const double orientations[4][4] = {
            { 1,  0, -1,  0}, // north vs. north, east, south, west
            { 0,  1,  0, -1}, // east vs. north, east, south, west
            {-1,  0,  1,  0}, // south vs. north, east, south, west
            { 0, -1,  0,  1}  // west vs. north, east, south, west
        };

    public:

        double directional_persistence_covariate(
            const CardinalDirection & x,
            const CardinalDirection & y
        ) {
            return orientations[static_cast<std::size_t>(x)][
                static_cast<std::size_t>(y)
            ];
        }

};

#endif