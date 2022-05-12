#include "spinodal_omp.h"
#include <gsl/gsl_rng.h>

// Generate array of random numbers between 0 and 1
void rand_ZeroToOne(int seed, double *random_ZeroToOne_array) {
        const gsl_rng_type *T;
        gsl_rng *r;

        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, seed);

        // Generate array of random numbers between 0 and 1
        for (int i = 0; i < num_points; i++) {
            random_ZeroToOne_array[i] = gsl_rng_uniform(r);
        }

        gsl_rng_free(r);

    }