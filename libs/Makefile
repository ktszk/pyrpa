#Makefile
#FC=ifort
FC=ifx
#FC=flang
#FC=gfortran
#FC=amdflang
ifeq ($(FC),ifort)
  LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core #-lmkl_avx
  FFLAGS= -O2 -mavx -qopenmp -IPO
else
  ifeq ($(FC),ifx)
    LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
    FFLAGS= -fPIC -O2 -mavx2 -fopenmp -finline-functions -fPIC
  else
    ifeq ($(FC),flang)
      LIBS= -lblis -lflame -lfftw3 -lm -lpthread
      FFLAGS= -O2 -mavx2 -fopenmp -finline-functions
    else
      #LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core
      LIBS= -lblis -lflame -lfftw3 -lm -lpthread
      FFLAGS= -O2 -mavx2 -fopenmp -finlne-functions
    endif
  endif
endif

FFLAGS2= -shared -fPIC
objmain= fmod.o fcond.o fchi.o fself.o
.SUFFIXES:
main: $(objmain)
	$(FC) $(FFLAGS) $(FFLAGS2) $(LIBS) -o fmod.so ${objmain}
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
.PHONY: clean
clean:
	rm -rf *.mod *.o *.so __pycache__
