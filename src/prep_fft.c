//
// Created by Jitin Nair on 24/01/20.
//
#include "spinodal_omp.h"

void prep_fft(double *k2, double *k4){

    double *kx, *ky;
    kx=(double*)malloc(sizeof(double)*num_points);
    ky=(double*)malloc(sizeof(double)*num_points);

    double delNx, delNy;
    delNx=2*M_PI/(Nx*dx);
    delNy=2*M_PI/(Ny*dy);

    //Periodic boundary conditions
    for(int i=0; i<Nx; i++){
        if(i < Nx/2)
            kx[i]=(double)(i)*delNx;
        else
            kx[i]=(double)(i-Nx)*delNx;
    }

    for(int j=0; j<Ny; j++){
        if(j<Ny/2)
            ky[j]=(double)(j)*delNy;
        else
            ky[j]=(double)(j-Ny)*delNy;
    }

    int ii;
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
                ii = j + i*Ny;
                k2[ii] = (kx[i] * kx[i]) + (ky[j] * ky[j]);
                k4[ii] = k2[ii] * k2[ii];
            }
        }

    free(kx);
    free(ky);
}
