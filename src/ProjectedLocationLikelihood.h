#ifndef MOVECON_PROJECTED_LOCATION_LIKELIHOOD_H
#define MOVECON_PROJECTED_LOCATION_LIKELIHOOD_H

#include <Rcpp.h>

/**
 * Implements "Materials and Methods" from McClintock et al. (2015)
 * doi: 10.1111/2041-210X.12311
*/
struct ProjectedLocationLikelihood {

    private:

        // bivariate normal parameters and constants
        double mu_easting, mu_northing, sd_easting, sd_northing, rho, rhosq_c, 
            lcst, conditional_scaling, conditional_sd;

    public:

        ProjectedLocationLikelihood(
            double easting, double northing, double semi_major, 
            double semi_minor, double orientation
        ) : mu_easting(easting), mu_northing(northing) {

            //
            // parameterize bivariate normal distribution
            //

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

            // log-normalizing constant for density evaluation
            lcst = - std::log(2 * M_PI * sd_easting * sd_northing)
                - 0.5 * std::log(rhosq_c); 

            //
            // parameterize constants that facilitate sampling
            //

            conditional_scaling = sd_northing / sd_easting * rho;
            conditional_sd = std::sqrt(rhosq_c) * sd_northing;
        }

        /**
         * Evaluate log-likelihood for a state
        */
       template<typename State>
       double dstate(const State & state) {
            // compute projected distances; scale wrt. uncertainty
            double zx = (
                state.properties.location->easting - mu_easting
            ) / sd_easting;
            double zy = (
                state.properties.location->northing - mu_northing
            ) / sd_northing;

            // sign the distances
            if(mu_easting < state.properties.location->easting) {
                zx *= -1;
            }
            if(mu_northing < state.properties.location->northing) {
                zy *= -1;
            }

            // quadratic form
            double q = zx * zx - 2 * rho * zx * zy + zy * zy;

            return - q / 2 / rhosq_c + lcst;
       }

        /**
         * Evaluate log-likelihood for a particle
        */
        template<typename Particle>
        double dparticle(const Particle & particle) {
            return dstate(*particle.state);
        }

        /**
         * Draw a sample from the parameterized distribution, using output 
         * parameters to write directly to pre-allocated memory
        */
        void sample(double & easting, double & northing) {
            // sample the joint distribution using normal conditional properties
            easting = R::rnorm(mu_easting, sd_easting);
            northing = R::rnorm(
                mu_northing + conditional_scaling * (easting - mu_easting),
                conditional_sd
            );
        }

};

/**
 * Create a family of location observation distributions from input vectors
*/
std::vector<ProjectedLocationLikelihood> LocationDistributionFamily(
    std::vector<double> eastings, std::vector<double> northings, 
    std::vector<double> semi_majors, std::vector<double> semi_minors,
    std::vector<double> orientations
);

#endif
