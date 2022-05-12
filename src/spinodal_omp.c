#include "spinodal_omp.h"

int Nx, Ny, Nz, num_points;
double dx, dy, dz, dt;

int main(int argc, char const *argv[]) {

    int sflag=0;            //enable serial execution by setting sflag=1

    //check if parallel execution enabled
    if(sflag==1) {
        omp_set_num_threads(1);
    }

    //start clock
    double start_time = omp_get_wtime();

    //Parallel FFTW Setup
    int status=fftw_init_threads();
    if (status!=0) {
        printf("Running in parallel with %d threads\n", omp_get_max_threads());
    }
    fftw_plan_with_nthreads(omp_get_max_threads());

    //read arguments and assign
    options_t options= {argv[1], argv[2]};

    //write input to file
    write_input_to_file(options);

    //declarations
    Nx=1024, Ny=Nx, Nz=0;
    dx=1.0, dy=1.0, dz=1.0;
    dt=0.5;

    double *k2, *k4, *conc_print;

    fftw_complex *conc, *conc_tilde, *free_energy, *free_energy_tilde;
    fftw_plan planForward, planBackward;

    //evolution constants
    double conc0=strtod(options.initial_conc, NULL),
            mobility=1.0,
            kappa=1.0,
            A=1.0;

    //output save constants
    int nstep=5000,
    iprint=500,
    istep=0;

    num_points = Nx * Ny;

    //output label
    char* label = malloc( sizeof(char) * (strlen(options.label)+1) );
    strcpy( label, options.label);

    //type of microstructure (iflag=1 is one grain setup for benchmarking)
    int iflag=2;

    //FFTW allocations
    conc=(fftw_complex*)fftw_malloc(num_points * sizeof(fftw_complex));
    conc_tilde=(fftw_complex*)fftw_malloc(num_points * sizeof(fftw_complex));

    free_energy=(fftw_complex*)fftw_malloc(num_points * sizeof(fftw_complex));
    free_energy_tilde= (fftw_complex*)fftw_malloc(num_points * sizeof(fftw_complex));

    //declare k-vectors for Fourier space
    k2=(double*)malloc(sizeof(double)*num_points);
    k4=(double*)malloc(sizeof(double)*num_points);

    //For RNG and write_to_VTK
    conc_print=(double*)malloc(sizeof(double) * num_points);

#pragma omp parallel default(none) shared(Nx, Ny, Nz, conc, conc_tilde,\
planForward, planBackward) private(istep)

    //FFTW plans
#pragma omp single nowait
    planForward=fftw_plan_dft_2d(Nx, Ny, conc, conc_tilde, -1, FFTW_MEASURE);
    planBackward=fftw_plan_dft_2d(Nx, Ny, conc_tilde, conc, +1, FFTW_MEASURE);

    //TODO: Check FFTW documentation. Should a 1D plan across total number of points be used here?

    //prepare microstructure
    prep_microstructure(iflag, conc, conc_print, conc0);

    //write initial concentration to file
#pragma omp single nowait
    write_to_VTK(istep , conc_print, label);

    //print completion status
    printf("Timestep %d completed\n", istep );

    //prepare kx and ky
    prep_fft(k2, k4);

    //time loop
    for(istep=1; istep<=nstep; istep++) {

        //calculate derivative of free_energy
#pragma omp for nowait
        for (int ii = 0; ii < num_points; ii++) {
                free_energy[ii] = 2 * A * creal(conc[ii]) * (1 - creal(conc[ii])) * (1 - 2 * creal(conc[ii]));
            }

        fftw_execute_dft(planForward, free_energy, free_energy_tilde);   //calculating free_energy_tilde
        fftw_execute_dft(planForward, conc, conc_tilde);    //calculating conc_tilde

        // calculation in fourier space
#pragma omp for nowait
        for (int ii = 0; ii < num_points; ii++) {
            double numer = dt * mobility * k2[ii] * free_energy_tilde[ii];
            double denom = 1.0 + dt * A * mobility * kappa * k4[ii];
            conc_tilde[ii] = (conc_tilde[ii] - numer) / denom + _Complex_I*0.0;
        }

        //coming to real space
        fftw_execute(planBackward);

        //time interval for saving the data
        if (istep % iprint == 0) {
            //write real part to conc[]
#pragma omp for nowait
            for (int ii = 0; ii < num_points; ii++) {
                conc_print[ii] = creal(conc[ii]) / (double) (num_points);
            }

            // write concentration to file
            write_to_VTK(istep, conc_print, label);

            //print completion status
            printf("Timestep %d completed\n", istep);
        }

        //exchanging the values
#pragma omp for nowait
        for (int ii = 0; ii < num_points; ii++) {
            conc[ii]= creal(conc[ii]) / (double)(num_points);
        }
    }

    //FFTW Free, Destroy and Cleanup
    free(label);
    free(k2);
    free(k4);
    free(conc_print);
    fftw_free(conc);
    fftw_free(conc_tilde);
    fftw_free(free_energy);
    fftw_free(free_energy_tilde);
    fftw_destroy_plan(planForward);
    fftw_destroy_plan(planBackward);
    fftw_cleanup();

    // end clock
    double time = omp_get_wtime() - start_time;
    printf("Elapsed: %f seconds\n", time);

    return 0;
}