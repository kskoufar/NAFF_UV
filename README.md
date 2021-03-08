# NAFF_UV
The NAFF_UV code is an implementation of the Jacques Laskar's NAFF algorithm.

This new version (NAFF_UV) can be seen as an updated and enhanced version of Jacques Laskar's original NAFF that was written in Fortran 77 and of a later version used in Bmad (https://github.com/MichaelEhrlichman/FortNAFF) and written in Fortran 90.

NAFF_UV is written in Fortran 2008 and there is a dependence on FFTW libraries (http://www.fftw.org/).

### For compilation use:
gfortran -o3 -I./ -o abc kind_accuracy_module.f08 numbers_constants_module.f08 brent_mnbrak_module.f08 naff_uv.f08 -L./ -lfftw3 -lm
 
### In lxplus use:
 lxplus7
### or: 
 scl enable devtoolset-7 bash

### References:
k. Skoufaris, "Non-linear dynamics modeling in accelerators with the use of symplectic integrators".

J. Laskar, "The chaotic motion of the solar system: A numerical estimate of the size of the chaotic zones".

J. Laskar, C. Froeschlé, and A. Celletti, "The measure of chaos by the numerical analysis of the fundamental frequencies. Application to the standard mapping".

Y. Papaphilippou, "Detecting chaos in particle accelerators through the frequency map analysis method".
