# pyrpa
 pyrpa is python script and fortran library for model calculation.
 pyrpa can calculate band structure, Dos, Fermi surface, conductivity, and etc. using tight binding model Hamiltonian.

# Requirements
 - fortran compiler (we check gcc 11.4 and newest aocc and ifx)
 - lapack and blas (we check OpenBLAS, MKL, AOCL(BLIS and Flame))
 - python3
 - python standard library
  -- __future__
  -- ctypes
 - python external packages
  -- numpy
  -- scipy
  -- scikit-image
  -- matplotlib

# Composition
 - main.py: main program written python
 - fmod.so: fortran library
 - plibs.py: python functions
 - flibs.py: wrapper functions of fmod.so for python
fmod.so,plibs.py and flibs.py are placed with the fortran codes in the libs dir.

# How to compile
 - clone this repository
 - enter the libs directory
 - edit Makefile, select compiler and math library.
 - do Make command then fmod.so is generated.