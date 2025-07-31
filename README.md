# pyrpa
pyrpa is python script and fortran library for model calculation.
pyrpa can calculate band structure, Dos, Fermi surface, conductivity, and etc. using tight binding model Hamiltonian.

# Requirements
- fortran compiler (we check gcc 13.3 and newest AMD Flang(ROCm and aocc) and ifx)
- lapack, blas, and FFTW (we check OpenBLAS, MKL, and AOCL(BLIS and Flame))
- python3 (we check 3.13.5)
- python standard library
  - ctypes
- python external packages
  - numpy
  - scipy
  - scikit-image
  - matplotlib

# Composition
- main.py: main program written python
- fmod.so: fortran library
- plibs.py: python functions
- flibs.py: wrapper functions of fmod.so for python

fmod.so, plibs.py and flibs.py are placed with the fortran codes in the libs dir.

# How to compile
- clone this repository
- enter the libs directory
- edit Makefile, select compiler and math library.
- do Make command then fmod.so is generated.
