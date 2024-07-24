#Makefile
#FC=ifort
FC=ifx
#FC=flang
#FC=amdflang
#FC=gfortran
ifeq ($(FC),ifort)
  LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core #-lmkl_avx
  FFLAGS= -O2 -mavx2 -qopenmp -IPO -fPIC -heap-arrays
  FFLAGS2= -shared
else
  ifeq ($(FC),ifx)
    LIBS= -qmkl
    FFLAGS= -O2 -mavx2 -fopenmp -finline-functions -fPIC -heap-arrays
    FFLAGS2= -shared
  else
    ifeq ($(FC),flang)
      LIBS= -lblis -lflame -lfftw3 -lm -lpthread
      FFLAGS= -O2 -mavx2 -fopenmp -finline-functions
      FFLAGS2= -shared -fPIC
    else
      ifeq ($(FC),amdflang)
        LIBS= -lblis -lflame -lfftw3 -lm -lpthread
        FFLAGS= -O2 -mavx2 -fopenmp -finline-functions -fPIC
        FFLAGS2= -shared
      else
        LIBS= -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -liomp5
        #LIBS= -lblis -lflame -lfftw3 -lm -lpthread
        #LIBS= -lopenblas -lfftw3
        FFLAGS= -O2 -mavx2 -fopenmp -finline-functions
        FFLAGS2= -shared -fPIC
      endif
    endif
  endif
endif

objmain= fmod.o fcond.o fchi.o fphi.o fself.o fimp.o fpade.o
.SUFFIXES:
main: $(objmain)
	$(FC) $(FFLAGS) $(FFLAGS2) $(LIBS) -o fmod.so ${objmain}
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
.PHONY: clean
clean:
	rm -rf *.mod *.o *.so __pycache__
