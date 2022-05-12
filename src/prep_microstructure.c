#include "spinodal_omp.h"

void prep_microstructure(int iflag, fftw_complex conc[num_points], double conc_print[num_points], double c0) {

  // random orientation
  if (iflag == 2) {
    double noise = 0.02;
    //get array of random numbers between 0 and 1 for setting initial microstructure
    double *random_ZeroToOne_array;
    random_ZeroToOne_array = (double *) malloc(sizeof(double) * num_points);
    rand_ZeroToOne(0, random_ZeroToOne_array);

#pragma omp for nowait
    for (int ii = 0; ii < num_points; ii++) {
      conc[ii] = c0 + noise * (0.5 - random_ZeroToOne_array[ii]);
      conc_print[ii] = creal(conc[ii]);
    }

    free(random_ZeroToOne_array);
  }
}
