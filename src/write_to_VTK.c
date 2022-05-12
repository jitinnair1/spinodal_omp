#include "spinodal_omp.h"

void write_to_VTK(int iprint, double conc[num_points], char label[30]) {

  // create output directory if not created
  char path[60]="../output/";
  struct stat st = {0};
  if (stat(path, &st) == -1) {
    mkdir(path, 0700);
  }

  strcat(path, label);
  strcat(path,"/");
  struct stat st2 = {0};
  if (stat(path, &st2) == -1) {
    mkdir(path, 0700);
  }

  // create folder within output directory named $label
  char time[30];
  sprintf(time, "time_%05d", iprint);
  char filename[strlen(path) + strlen(time) + 5];
  strcpy(filename, path);
  strcat(filename, time);
  strcat(filename, ".vtk");

  FILE *fp;
  fp = fopen(filename,"w");

  // write header of VTK file
  fprintf(fp, "# vtk DataFile Version 2.0\n");
  fprintf(fp, "time_10.vtk\n");
  fprintf(fp, "ASCII\n");

  // co-ordinates of grid points
  fprintf(fp, "DATASET STRUCTURED_POINTS\n");
  int tempNz = 1;
  fprintf(fp, "DIMENSIONS %6d %6d %6d\n", Nx, Ny, tempNz);
  fprintf(fp, "ORIGIN 0 0 0\n");
  fprintf(fp, "SPACING %f %f %f\n", dx, dy, dz );

  // grid point values for C
  fprintf(fp,"POINT_DATA %6d\n", num_points);
  fprintf(fp,"SCALARS CONC  float  1\n");
  fprintf(fp,"LOOKUP_TABLE default\n");
    for (int j = 0; j < Ny; j++) {
      for (int i = 0; i < Nx; i++) {
        int ii = j + i * Ny;
        fprintf(fp, "%lf\n", conc[ii]);
      }
    }
  fclose(fp);
}
