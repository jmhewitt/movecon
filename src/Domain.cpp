#include "Domain.h"

RookDirectionalStatespace::RookDirectionalStatespace(
    const Rcpp::NumericVector & eastings, 
    const Rcpp::NumericVector & northings,
    Rcpp::NumericMatrix & covariates,
    Rcpp::NumericVector & linear_constraint
) {

    /* 
      Spatial grid format reference: 
      https://cran.r-project.org/web/packages/stars/vignettes/stars4.html 
    */

    // direction of north/east with respect to grid
    int north_step = northings[1] - northings[0] > 0 ? 1 : -1;
    int east_step = eastings[1] - eastings[0] > 0 ? 1 : -1;

    // extract grid metadata
    std::size_t eastings_len = eastings.size();
    std::size_t northings_len = northings.size();
    double min_easting, max_easting, min_northing, max_northing;
    if(north_step == 1) {
        min_northing = northings[0];
        max_northing = northings[northings_len-1];
    } else {
        max_northing = northings[0];
        min_northing = northings[northings_len-1];
    }
    if(east_step == 1) {
        min_easting = eastings[0];
        max_easting = eastings[eastings_len-1];
    } else {
        max_easting = eastings[0];
        min_easting = eastings[eastings_len-1];
    }
    
    // prepare to associate covariate information with locations
    std::size_t p = covariates.rows();
    auto covariates_it = covariates.begin() - p;

    // access the linear constraint as a vector
    Eigen::Map<Eigen::VectorXd> linear_constraint_vec(
        linear_constraint.begin(), linear_constraint.size()
    );

    // build grid
    for(std::size_t j = 0; j < northings_len; ++j) {
        for(std::size_t i = 0;  i < eastings_len; ++i) {

            // update pointers
            covariates_it += p;

            // skip location if linear constraint is not satisfied
            Eigen::Map<Eigen::VectorXd> x(covariates_it, p);
            if(linear_constraint_vec.dot(x) < 0)
                continue;
        
            // initialize grid cell; address never changes b/c grid is std::map
            Location & cell = grid[LocationIndices(i,j)];

            // transfer information to grid cell
            cell.easting = eastings[i];
            cell.northing = northings[j];
            new (&(cell.x)) Eigen::Map<Eigen::VectorXd>(covariates_it, p);
        }
    }

    // initialize states associated with grid cells
    for(auto& map_entry : grid) {

        std::size_t easting_ind = map_entry.first.first;
        std::size_t northing_ind = map_entry.first.second;
        Location & cell = map_entry.second;

        // define neighbor keys
        LocationIndices southern_neighbor(
            easting_ind, northing_ind - north_step
        );
        LocationIndices northern_neighbor(
            easting_ind, northing_ind + north_step
        );
        LocationIndices western_neighbor(
            easting_ind - east_step, northing_ind
        );
        LocationIndices eastern_neighbor(
            easting_ind + east_step, northing_ind
        );
        
        // cell can be reached from a more southern location
        if(grid.count(southern_neighbor) > 0) {
            StateType state;
            state.properties.location = &cell;
            state.properties.last_movement_direction = north;
            states[StateKey(north,easting_ind,northing_ind)] = state;
        }

        // cell can be reached from a more northern location
        if(grid.count(northern_neighbor) > 0) {
            StateType state;
            state.properties.location = &cell;
            state.properties.last_movement_direction = south;
            states[StateKey(south,easting_ind,northing_ind)] = state;
        }

        // cell can be reached from a more western location
        if(grid.count(western_neighbor) > 0) {
            StateType state;
            state.properties.location = &cell;
            state.properties.last_movement_direction = east;
            states[StateKey(east,easting_ind,northing_ind)] = state;
        }

        // cell can be reached from a more eastern location
        if(grid.count(eastern_neighbor) > 0) {
            StateType state;
            state.properties.location = &cell;
            state.properties.last_movement_direction = west;
            states[StateKey(west,easting_ind,northing_ind)] = state;
        }
    }

    // initialize connections between states 
    // (i.e., link allowable movement combinations on the grid)
    for(auto& map_entry : states) {

        std::size_t easting_ind = std::get<1>(map_entry.first);
        std::size_t northing_ind = std::get<2>(map_entry.first);
        StateType & state = map_entry.second;

        // define movement keys
        StateKey eastern_movement = StateKey(
            east, easting_ind + east_step, northing_ind
        );
        StateKey western_movement = StateKey(
            west, easting_ind - east_step, northing_ind
        );
        StateKey northern_movement = StateKey(
            north, easting_ind, northing_ind + north_step
        );
        StateKey southern_movement = StateKey(
            south, easting_ind, northing_ind - north_step
        );

        // allow movement to the east
        if(states.count(eastern_movement) > 0) {
            // state associated with movement to eastern neighbor
            StateType & nbr = states.at(eastern_movement);
            // forward link
            state.to.insert(&nbr);
            // backward link
            nbr.from.insert(&state);
        }

        // allow movement to the west
        if(states.count(western_movement) > 0) {
            // state associated with movement to western neighbor
            StateType & nbr = states.at(western_movement);
            // forward link
            state.to.insert(&nbr);
            // backward link
            nbr.from.insert(&state);
        }

        // allow movement to the south
        if(states.count(southern_movement) > 0) {
            // state associated with movement to southern neighbor
            StateType & nbr = states.at(southern_movement);
            // forward link
            state.to.insert(&nbr);
            // backward link
            nbr.from.insert(&state);
        }

        // allow movement to the north
        if(states.count(northern_movement) > 0) {
            // state associated with movement to northern neighbor
            StateType & nbr = states.at(northern_movement);
            // forward link
            state.to.insert(&nbr);
            // backward link
            nbr.from.insert(&state);
        }
    }
}

/**
 * Create a linked-list representation of a discrete state space for persistent
 * movement with rook adjacencies.  Returns an Rcpp::XPtr to the linked-list in
 * C++.
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
// [[Rcpp::export]]
Rcpp::XPtr<RookDirectionalStatespace> build_statespace(
    Rcpp::NumericVector & eastings, 
    Rcpp::NumericVector & northings, 
    Rcpp::NumericMatrix & covariates,
    Rcpp::NumericVector & linear_constraint
) {
    RookDirectionalStatespace * statespace = new RookDirectionalStatespace(
        eastings, northings, covariates, linear_constraint
    );
    return Rcpp::XPtr<RookDirectionalStatespace>(statespace, true);
}

/**
 * Format a location object for viewing within R
*/
Rcpp::List format_location(const Location & location) {
    return Rcpp::List::create(
        Rcpp::Named("easting") = location.easting, 
        Rcpp::Named("northing") = location.northing,
        Rcpp::Named("covariates") = location.x
    );
}

/**
 * View a location from a CTDS domain object
 * 
 * @param statespace Object constructed from \code{build_statespace}
 * @param easting_ind 0-based index for location's easting coordinate
 * @param northing_ind 0-based index for location's northing coordinate
*/
// [[Rcpp::export]]
Rcpp::List extract_statespace_location(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    std::size_t easting_ind, 
    std::size_t northing_ind
) {
    Location & location = statespace->grid.at(
        RookDirectionalStatespace::LocationIndices(easting_ind, northing_ind)
    );
    return format_location(location);
}

/**
 * Convert a direction string to a CTDS domain object type
*/
RookDirectionalStatespace::Direction stringToDirection(
    const std::string & direction
) {
    if(direction.compare("north") == 0) {
        return RookDirectionalStatespace::Direction::north;
    } else if(direction.compare("east") == 0) {
        return RookDirectionalStatespace::Direction::east;
    } else if(direction.compare("south") == 0) {
        return RookDirectionalStatespace::Direction::south;
    } else if(direction.compare("west") == 0) {
        return RookDirectionalStatespace::Direction::west;
    }
    Rcpp::stop("Argument direction contains an invalid value");
}

/**
 * Format a CTDS domain object direction type for printing
*/
std::string directionToString(
    const RookDirectionalStatespace::Direction & direction
) {
    switch(direction) {
        case RookDirectionalStatespace::Direction::north:
            return "north";
        case RookDirectionalStatespace::Direction::east:
            return "east";
        case RookDirectionalStatespace::Direction::south:
            return "south";
        case RookDirectionalStatespace::Direction::west:
            return "west";
    }
    Rcpp::stop("String conversion is not defined for argument direction");
}

/**
 * Format the location-direction pair information from a CTDS domain state 
 * object for viewing in R
*/
Rcpp::List format_directional_persistence(
    const RookDirectionalStatespace::StateType & state
) {
    return Rcpp::List::create(
        Rcpp::Named("last_movement_direction") = directionToString(
            state.properties.last_movement_direction
        ),
        Rcpp::Named("location") = format_location(*state.properties.location)
    );
}

/**
 * Format a CTDS domain state object for viewing in R
*/
Rcpp::List format_state(const RookDirectionalStatespace::StateType & state) {
    Rcpp::List res = format_directional_persistence(state);
    Rcpp::List to = Rcpp::List::create();
    for(auto destination : state.to)
        to.push_back(format_directional_persistence(*destination));
    res["to"] = to;
    Rcpp::List from = Rcpp::List::create();
    for(auto source : state.from)
        from.push_back(format_directional_persistence(*source));
    res["from"] = from;
    return res;
}

/**
 * View a state from a CTDS domain object
 * 
 * @param statespace Object constructed from \code{build_statespace}
 * @param last_movement_direction string "north", "east", "south", or "west" to
 *  identify the state to view
 * @param easting_ind 0-based index for location's easting coordinate
 * @param northing_ind 0-based index for location's northing coordinate
*/
// [[Rcpp::export]]
Rcpp::List extract_statespace_state(
    Rcpp::XPtr<RookDirectionalStatespace> statespace, 
    std::string last_movement_direction,
    std::size_t easting_ind, 
    std::size_t northing_ind
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
    return format_state(state);   
}
