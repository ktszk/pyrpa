#Makefile
#FC=ifort
#FC=ifx
FC=flang
#FC=gfortran
#FC=amdflang
ifeq ($(FC),ifort)
  LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core #-lmkl_avx
  parallel= -mavx -qopenmp
else
  ifeq ($(FC),ifx)
    LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
    parallel= -mavx2 -fopenmp
  else
    ifeq ($(FC),flang)
      LIBS= -lblis -lflame -lm -lpthread
      parallel= -mavx2 -fopenmp
    else
      LIBS= -lblis -lflame -lm -lpthread
      parallel= -mavx2 -fopenmp
    endif
  endif
endif
.SUFFIXES:
main:
	$(FC) -O2 $(parallel) -shared -fPIC $(LIBS) -o fmod.so fmod.f90
