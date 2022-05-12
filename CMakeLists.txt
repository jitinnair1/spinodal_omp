cmake_minimum_required(VERSION 3.14)
project(spinodal_omp)

#Set C Standard
set(CMAKE_C_STANDARD 99)

#Change binary output folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ../bin)

# Define flags for specific builds
set(CMAKE_C_FLAGS_DEBUG "-Wall -g -fopenmp -O2")

#NOTE: Debug should be run using the -O0 else values of variables are removed

# GCC include directory on macOS which includes openmp.h, openacc.h etc.
include_directories(/usr/local/Cellar/gcc)
include_directories(/usr/local/include)

# Link Libraries
link_directories(/usr/local/lib/)

# Source Files
include_directories(src)

add_executable(${PROJECT_NAME}
        src/prep_microstructure.c
        src/rand_ZeroToOne.c
        src/spinodal_omp.c
        src/spinodal_omp.h
        src/write_to_VTK.c
        src/write_init_conc.c
        src/prep_fft.c
        src/write_input_to_file.c)

# Link Library to Target
target_link_libraries(${PROJECT_NAME} fftw3_omp fftw3 m gsl gslcblas)





