#ifndef MOVECON_PARTICLE_FILTER_H
#define MOVECON_PARTICLE_FILTER_H

#include <Rcpp.h>

#include "log_add.h"

/**
 * The BootstrapParticleFilter uses an observer concept to export filtering 
 * distributions and incremental likelihood contributions.  The NullObsrever 
 * fulfills the observer requirement for BootstrapParticleFilter objects when 
 * the filtering distributions do not need to be exported.
*/
template<typename Particle>
struct NullObserver {
    void operator()(const std::vector<Particle> & particles, double ll) { }
};

/**
 * Store filtering distributions in an array
*/
template<typename Particle> 
struct FilterObserver {
    std::vector<std::vector<Particle>> particle_distributions;
    void operator()(const std::vector<Particle> & particles, double ll) { 
        particle_distributions.emplace_back(particles);
    }
};

template<
    typename Particle, 
    typename ProposalDistributionSequence, 
    typename LikelihoodSequence,
    typename Observer = NullObserver<Particle>
> 
class BootstrapParticleFilter {

    private:

        std::vector<Particle> particles_init;

        // convert a pointer to a reference if needed
        template<typename T> 
        T& asReference(std::unique_ptr<T> & x) { return  *x; }
        template<typename T> 
        T& asReference(T * x) { return  *x; }
        template<typename T> 
        T& asReference(T & x) { return  x; }

    public:

        ProposalDistributionSequence * proposal_distributions;
        LikelihoodSequence * likelihoods;

        /**
         * @param particles initial particles
        */
        BootstrapParticleFilter(const std::vector<Particle> & particles) : 
            particles_init(particles) { }

        /**
         *  Particle filter approximation to marginal log-likelihood using 
         *  a new, default observer
        */
        double marginal_ll() {
            Observer o;
            return marginal_ll(o);
        }

        /**
         * Particle filter approximation to marginal log-likelihood, implemented
         * via Algorithm 1 (Bootstrap filter) of Michaud et. al. (2021, doi:
         * 10.18637/jss.v100.i03)
        */
        double marginal_ll(Observer & observer) {

            // initialize log-likelihood
            double ll = 0;

            // particle filter size
            std::size_t M = particles_init.size();
            double log_M = std::log(M);

            // set initial particle values (line 2)
            std::vector<Particle> particles_A = particles_init;
            std::vector<Particle>* active_particles = &particles_A;

            // compute initial weights (line 3)
            double log_uniform_weight = -std::log(M);

            // prepare container for unnormalized weights (line 8)
            std::vector<double> log_unnormalized_weights(M);

            // prepare container for resampling (line 14)
            std::vector<Particle> particles_B;
            particles_B.reserve(M);
            std::vector<Particle>* resampled_particles = &particles_B;

            // iterate over observations (line 5)
            auto proposal_distn = proposal_distributions->begin();
            auto proposal_distn_end = proposal_distributions->end();
            auto likelihood = likelihoods->begin();
            for(; proposal_distn != proposal_distn_end; ++proposal_distn) {

                // evaluate proposal distributions and importance weights
                auto particle = active_particles->begin();
                auto particle_end = active_particles->end();
                auto log_w = log_unnormalized_weights.begin();
                for(; particle != particle_end; ++particle) {
                    // sample from proposal distribution (line 7)
                    proposal_distn->propose(*particle);
                    // compute log-importance weight (line 8)
                    // Note: weight will always be uniform
                    *(log_w++) = asReference(*likelihood).dparticle(*particle) + 
                        log_uniform_weight;
                }

                // normalize resampling weights and resample (lines 11, 14, 15)
                std::size_t nresample = M;
                particle = active_particles->begin();
                resampled_particles->clear();
                double log_mass = log_sum(log_unnormalized_weights);
                log_w = log_unnormalized_weights.begin();
                auto log_w_stop = log_unnormalized_weights.end() - 1;
                for(; log_w != log_w_stop; ++log_w) {
                    // normalized resampling weight (line 11)
                    double p = std::exp(*log_w - log_mass);
                    // resample particles using R::multinom strategy
                    // (see R source code: R-XXX/src/nmath/rmultinom.c)
                    if(p != 0) {
                        // number of times to use particle (line 14)
                        std::size_t n = static_cast<std::size_t>(
                            R::rbinom(nresample, p)
                        );
                        // transfer n copies of current particle (line 15)
                        if(n > 0) {
                            nresample -= n;
                            resampled_particles->insert(
                                resampled_particles->end(), n, *particle
                            );
                        }
                    }
                    // iterate particle
                    ++particle;
                }
                
                // transfer final particle, if needed
                if(nresample > 0) {
                    resampled_particles->insert(
                        resampled_particles->end(), nresample, *particle
                    );
                }

                // aggregate likelihood mass (line 18)
                double ll_t = log_mass - log_M;
                ll += ll_t;

                // update particles
                std::swap(active_particles, resampled_particles);

                // provide opportunity to export filtering distributions, etc.
                observer(*active_particles, ll_t);

                // increment likelihood
                ++likelihood;

            } // iterate over observations (line 5)

            return ll;
        } // marginal_ll()

};

#endif
