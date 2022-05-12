#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gperftools/profiler.h>
#include <omp.h>

typedef struct {
    char const*          initial_conc;
    char const*          label;
} options_t;

extern int Nx, Ny, Nz, num_points;
extern double dx, dy, dz, dt;

void prep_microstructure(int iflag, fftw_complex conc[num_points], double conc_print[num_points], double c0) ;

void rand_ZeroToOne(int seed, double *random_ZeroToOne_array) ;

void write_to_VTK(int iprint, double conc[num_points], char label[30]);

void prep_fft(double *k2, double *k4);

void write_input_to_file(options_t options);

