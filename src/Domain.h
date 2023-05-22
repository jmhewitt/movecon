/**
 * Objects that define the state space for a CTDS model with directional 
 * persistence
*/

#ifndef MOVECON_DOMAIN_H
#define MOVECON_DOMAIN_H

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

/**
 * Basic description for spatial information
*/
struct Location {

    // coordinates
    double easting, northing;

    // covariates
    Eigen::Map<Eigen::VectorXd> x{nullptr, 0};

};

/**
 * Basic description for directional persistence, associating a location with 
 * the most recent direction of travel taken to get to the location, e.g., to 
 * a location from the north.  A DirectionalPersistence object can serve as a 
 * state in a movement model with directional persistence.
*/
template<typename Direction, typename Location>
struct DirectionalPersistence {

    Direction last_movement_direction;

    Location location;

};

/**
 * State objects will not be 1:1 with Location objects for models that include 
 * directional persistence.  For such models, State objects will encode the 
 * current location as well as the most recent direction of travel.
*/
template<typename Properties>
struct State {

    typedef State<Properties> SelfType;

    // application-specific information about the state's interpretation
    Properties properties;

    // states to which transitions can be made, e.g., to the north
    std::set<SelfType*> to;

    // states from which transitions originate, e.g., from the north
    std::set<SelfType*> from;
};

/**
 * Assumes regular grid without affine transformations
*/
struct RookDirectionalStatespace {
    
    // define neighborhoods and movement with respect to rook moves
    enum Direction{north, east, south, west};
    typedef DirectionalPersistence<Direction, Location*> LocalMovement;
    
    // grid is a collection of locations
    typedef std::pair<std::size_t, std::size_t> LocationIndices;
    std::map<LocationIndices, Location> grid;

    // statespace is collection of possible transitions between grid cells
    typedef std::tuple<Direction, std::size_t, std::size_t> StateKey;
    typedef State<LocalMovement> StateType;
    std::map<StateKey, StateType> states;

    /**
     * Linked-list representation of a discrete state space for persistent 
     * movement with rook adjacencies.
     * 
     * Each state pairs a spatial location with the direction of movement used 
     * to arrive at that location (i.e., N,E,S,W for rook adjacencies).  Links
     * between states facilitate fast transitions on the state space.  Each 
     * location is associated with a covariate vector.
     * 
     * TODO: Add constraints to the spatial domain.
     * 
     * @param eastings vector of easting coordinates in monotonic order, either
     *   increasing or decreasing
     * @param northings vector of northing coordinates in monotonic order, either 
     *   increasing or decreasing
     * @param covariates matrix in column-major storage format in which each 
     *  column defines the covariates for a spatial location. A nested loop 
     *  of the northing and easting vectors forms the spatial ordering of the 
     *  columns. The northing coordinates are the outer loop, and the easting
     *  coordinates are the inner loop.  The first block of columns iterates
     *  over all easting coordinates for the first northing coordinate, etc.
    */
    RookDirectionalStatespace(
        const Rcpp::NumericVector & eastings,
        const Rcpp::NumericVector & northings,
        Rcpp::NumericMatrix & covariates
    );
};

#endif //MOVECON_DOMAIN_H
