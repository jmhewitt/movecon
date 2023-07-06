#ifndef MOVECON_PROJECTED_LOCATION_LIKELIHOOD_H
#define MOVECON_PROJECTED_LOCATION_LIKELIHOOD_H

#include <Rcpp.h>

/**
 * Implements "Materials and Methods" from McClintock et al. (2015)
 * doi: 10.1111/2041-210X.12311
*/
template<typename Particle>
struct ProjectedLocationLikelihood {

    private:

        // bivariate normal parameters and constants
        double mu_easting, mu_northing, sd_easting, sd_northing, rho, rhosq_c, 
            lcst;

    public:

        ProjectedLocationLikelihood(
            double easting, double northing, double semi_major, 
            double semi_minor, double orientation
        ) : mu_easting(easting), mu_northing(northing) {

            // parameterize bivariate normal distribution

            double Mtsq_half = semi_major * semi_major / 2;
            double mtsq_half = semi_minor * semi_minor / 2;

            double c = orientation / 180 * M_PI;
            double cos_c = std::cos(c);
            double cos2_c = cos_c * cos_c;
            double sin_c = std::sin(c);
            double sin2_c = sin_c * sin_c;

            sd_easting = std::sqrt(Mtsq_half * sin2_c + mtsq_half * cos2_c);
            sd_northing = std::sqrt(Mtsq_half * cos2_c + mtsq_half * sin2_c);
            rho = (Mtsq_half - mtsq_half) * cos_c * sin_c / sd_easting / 
                sd_northing;

            rhosq_c = 1 - rho * rho;

            // log-normalizing constant
            lcst = - std::log(2 * M_PI * sd_easting * sd_northing)
                - 0.5 * std::log(rhosq_c); 
        }

        double operator()(const Particle & particle) {

            // compute projected distances; scale wrt. uncertainty
            double zx = (
                particle.state->properties.location->easting - mu_easting
            ) / sd_easting;
            double zy = (
                particle.state->properties.location->northing - mu_northing
            ) / sd_northing;

            // sign the distances
            if(mu_easting < particle.state->properties.location->easting) {
                zx *= -1;
            }
            if(mu_northing < particle.state->properties.location->northing) {
                zy *= -1;
            }

            // quadratic form
            double q = zx * zx - 2 * rho * zx * zy + zy * zy;

            return - q / 2 / rhosq_c + lcst;
        }

};

/**
 * Create a family of location observation distributions from input vectors
*/
template<typename Particle>
std::vector<ProjectedLocationLikelihood<Particle>> LocationDistributionFamily(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations
) {
    std::vector<ProjectedLocationLikelihood<Particle>> family;
    family.reserve(eastings.size());
    
    auto eastings_it = eastings.begin();
    auto northings_it = northings.begin();
    auto semi_majors_it = semi_majors.begin();
    auto semi_minors_it = semi_minors.begin();
    auto orientations_it = orientations.begin();

    auto eastings_end = eastings.end();
    for(; eastings_it != eastings_end; ++eastings_it) {
        family.emplace_back(
            ProjectedLocationLikelihood<Particle>(
                *eastings_it, *(northings_it++), *(semi_majors_it++), 
                *(semi_minors_it++), *(orientations_it++)
            )
        );
    }

    return family;
}

#endif
