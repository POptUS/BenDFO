# BenDFO
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)    [![miss_hit](https://github.com/POptUS/BenDFO/actions/workflows/miss_hit.yml/badge.svg)](https://github.com/POptUS/BenDFO/actions/workflows/miss_hit.yml)

## Benchmarking Derivative-Free Optimization Algorithms 

This repository contains source code and updates from the codes originally shared at https://www.mcs.anl.gov/~more/dfo/.
Ways to contribute; extensions to other problems and languages; and additional resources will be updated here.

---

The codes in this repository are based on the supplemental information for the paper:
- [[1](#pap1)] Benchmarking Derivative-Free Optimization Algorithms by J.J. Moré and S.M. Wild. *SIAM J. Optimization*, Vol. 20 (1), pp. 172-191, 2009. doi:[10.1137/080724083](https://doi.org/10.1137/080724083)
   

### Benchmark Problems
The following source files were used to define the benchmark problems in [[1](https://github.com/POptUS/BenDFO#pap1)]. Fortran files are found in `fortran/` and Octave/Matlab files are found in `m/`.

- `calfun` [[Fortran](fortran/calfun.f)]  [[Matlab/Octave](m/calfun.m)]:
  Code for evaluating the 22 CUTEr problems considered (needs dfovec).
- `dfovec` [[Fortran](fortran/dfovec.f)]  [[Matlab/Octave](m/dfovec.m)]:
  Code producing component vectors for the 22 CUTEr problems considered.
- `dfoxs` [[Fortran](fortran/dfoxs.f)]  [[Matlab/Octave](m/dfoxs.m)]:
  Code specifying the standard starting points.
- `dfo.dat` [Data file](data/dfo.dat):
  Data file specifying the benchmark problem set P through the integer parameters (nprob, n, m, ns). 

### Plotting the Profiles
We provide the following Octave/Matlab files for producing basic data and performance profiles from data:

- `data_profile` [[Matlab/Octave](profiling/data_profile.m)]:
    Code for plotting a basic data profile.
- `perf_profile` [[Matlab/Octave](profiling/perf_profile.m)]:
    Code for plotting a basic performance profile.
    
Updates and other languages will be reflected here. Please also see:

- [Julia package for data and performance profiles](https://github.com/JuliaSmoothOptimizers/BenchmarkProfiles.jl)


### Derivatives and Testing
Autodiff versions (using [adimat](https://www.informatik.tu-darmstadt.de/sc/res/sw/adimat/index.en.jsp)) of the problem derivatives are also included.

A sample calling script for testing and to see these derivative capabilities is provided in `testbendfo.m` [[Matlab/Octave](profiling/testbendfo.m)].

### Sample Solvers

Many of the original solvers considered in [[1](https://github.com/POptUS/BenDFO#pap1)] have seen significant refinements. Links to the original solvers considered are at https://www.mcs.anl.gov/~more/dfo/shootout.html 

## Contributing to BenDFO

Contributions are welcome in a variety of forms; please see [CONTRIBUTING](CONTRIBUTING.rst).

## License 

All code included in BenDFO is open source, with the particular form of license contained in the top-level 
subdirectories.  If such a subdirectory does not contain a LICENSE file, then it is automatically licensed 
as described in the otherwise encompassing BenDFO [LICENSE](/LICENSE).  


## Resources

To seek support or report issues, e-mail:

 * ``poptus@mcs.anl.gov``
