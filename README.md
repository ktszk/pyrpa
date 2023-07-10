# pyrpa
pyrpa is python script and fortran library for model calculation.
pyrpa can calculate band structure, Dos, Fermi surface, conductivity, and etc. using tight binding Hamiltonian.

# Requirements
- fortran compiler (we check  gcc 11.3 and llvm 14.06)
- lapack and blas (we check OpenBLAS, MKL, AOCL(BLIS and Flame))
- python3
- scipy
- scikit-image
- matplotlib

# How to compile
- open Makefile and select compiler.
- do Make command then fmod.so is generated.
- run main.py 