/**
 * Wrap boost methods to search a statespace for state and location entries
*/

#ifndef MOVECON_DOMAINSEARCH_H
#define MOVECON_DOMAINSEARCH_H

#include <Rcpp.h>

// [[Rcpp::depends(BH)]]

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "Domain.h"

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

template<typename Statespace>
class StatespaceSearch {

    private:

        // data structure to facilitate spatial searches of locations
        typedef bg::model::point<double, 2, bg::cs::cartesian> point;
        typedef std::pair<point, Location*> rtree_value;
        bgi::rtree<rtree_value, bgi::rstar<16>> domain_location_tree;

    public:

        // reverse-lookup to index into states by their location
        typedef typename Statespace::StateType StateType;
        std::map<Location*, std::set<StateType*>> states_by_location;  

        StatespaceSearch(Statespace & statespace) {
            // populate R-tree and reverse lookup structures with domain info.
            auto state = statespace.states.begin();
            auto state_end = statespace.states.end();
            for(; state != state_end; ++state) {
                Location * location = state->second.properties.location;
                // add location to reverse lookup for states
                states_by_location[location].insert(&state->second);
                // add location to rtree
                point p(location->easting, location->northing);
                domain_location_tree.insert(std::make_pair(p, location));
            }
        }

        /**
         * Return the location closest to the given coordinates
        */
        Location* map_location(double easting, double northing) {
            std::vector<rtree_value> res;
            domain_location_tree.query(
                bgi::nearest(point(easting, northing), 1), 
                std::back_inserter(res)
            );
            return res[0].second;
        }

};

#endif
