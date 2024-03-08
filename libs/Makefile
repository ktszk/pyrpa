#Makefile
#FC=ifort
FC=ifx
#FC=flang
#FC=gfortran
#FC=amdflang
ifeq ($(FC),ifort)
  LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core #-lmkl_avx
  parallel= -mavx -qopenmp -IPO
else
  ifeq ($(FC),ifx)
    LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
    parallel= -mavx2 -fopenmp -finline-functions
  else
    ifeq ($(FC),flang)
      LIBS= -lblis -lflame -lfftw3 -lm -lpthread
      parallel= -mavx2 -fopenmp
    else
      #LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
      LIBS= -lblis -lflame -lfftw3 -lm -lpthread
      parallel= -mavx2 -fopenmp -finlne-functions
    endif
  endif
endif

FFLAGS= -O2 $(parallel) -shared -fPIC
objmain= fmod.o fcond.o fchi.o
.SUFFIXES:
main: $(objmain)
	$(FC) $(FFLAGS) $(LIBS) -o fmod.so ${objmain}
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
.PHONY: clean
clean:
	rm -rf *.mod *.so __pycache__
