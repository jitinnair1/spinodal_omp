# About

This repo has a simple implementation of the Cahn-Hilliard equation in C parallelised using OpenMP. I am posting my 
code here so that somebody starting out could use this as a template to build upon or to understand how the 
Cahn-Hilliard equation is solved using OpenMP parallelisation.

# Build

Use the following commands in your terminal:

```bash
git clone git@github.com:jitinnair1/hello-again-phasefield.git
cd hello-again-phasefield && mkdir build && cd build
cmake ..
make
```

# Usage

Once the project is built, you will find the executables under the bin folder. To run the executable use ./ followed 
by the name of the executable like:

```bash
./spinodal_omp X Y
```

where X and Y are the initial concentration and label for a given run. So if you want to set the initial concentration to be `0.5` and the use the label `Test`, you would run the 
executable as

```bash
./spinodal_omp 0.50 Test
```

This will generate a folder called `Test` within an output folder with .vtk files which can be opened 
using ParaView

# Issues

Please start a new discussion or issue if you encounter problems
