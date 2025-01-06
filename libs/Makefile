#Makefile
#FC=ifort
FC=ifx
#FC=flang
#FC=amdflang
#FC=gfortran
SL=MKL
#SL=AOCL
#SL=OpenBLAS
ifeq ($(FC),ifort)
  LIBS=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core #-lmkl_avx
  FFLAGS= -O2 -mavx2 -qopenmp -IPO -fPIC -heap-arrays
  FFLAGS2= -shared
else
  ifeq ($(FC),ifx)
    LIBS= -qmkl
    FFLAGS= -O2 -mavx2 -qopenmp -finline-functions -fPIC -heap-arrays -assume nan_compares
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
        ifeq ($(SL),MKL)
          LIBS= -lmkl_rt
        else
          ifeq ($(SL),AOCL)
            LIBS= -lblis -lflame -lfftw3 -lm -lpthread
          else
            LIBS= -lopenblas -lfftw3
          endif
        endif
        FFLAGS= -O2 -mavx2 -fopenmp -finline-functions
        FFLAGS2= -shared -fPIC
      endif
    endif
  endif
endif

srcdir= src
objdir= obj
flibso=fmod.so

objmain= $(addprefix $(objdir)/, fmod.o fcond.o fUJ.o fchi.o fphi.o fself.o fimp.o fpade.o feliash.o finv.o)
.SUFFIXES:
all: $(objdir) $(objmain) $(flibso)
$(flibso): $(objmain)
	$(FC) $(FFLAGS) $(FFLAGS2) $(LIBS) -o fmod.so $(objmain)
$(objdir):
	mkdir -p $(objdir)
$(objdir)/%.o: $(srcdir)/%.f90
	$(FC) $(FFLAGS) -o $@ -c $<
.PHONY: clean
clean:
	rm -rf *.mod *.f90 $(objdir) *.so __pycache__
