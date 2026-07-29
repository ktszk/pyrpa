#Makefile
# Compiler to use; override with `make FC=...` or environment variable.
# valid values: ifort, ifx, flang, amdflang, gfortran
#
# Default toolchain is ifx+MKL (the toolchain the shipped libfmod.so is built
# with). On these AMD CPUs `-march=<ARCH>` makes ifx emit vector code with NO
# Intel CPUID dispatch (verified), so it runs at full width on AMD.
#
# flang/amdflang (ROCm's LLVM flang) + AOCL is a fully supported alternative:
#   module load AMD/ROCm AMD/AOCL && make FC=flang
# AOCL supplies BLAS/LAPACK (libblis/libflame); the FFT comes from the bundled
# pocketfft (see the FFT block below for why AOCL's own FFTW3 is not used).
#
# --- Target CPU architecture -------------------------------------------------
# ARCH selects the -march target. Default znver3 is the lowest common
# denominator across the machines this runs on, so ONE binary works everywhere:
#   znver3  Zen 3  (e.g. Ryzen 7 5700G) -- AVX2, portable, DEFAULT
#   znver4  Zen 4  (e.g. Ryzen 9 7900)  -- AVX-512, faster but SIGILLs on Zen 3
# Build the portable default with `make`; on a Zen 4 box get full AVX-512 with:
#   make ARCH=znver4
ARCH ?= znver3
#
# For MKL BLAS/FFT to also use wider kernels on these AMD CPUs, run with e.g.
#   export MKL_ENABLE_INSTRUCTIONS=AVX512   (Zen 4)   or  AVX2 (Zen 3)
#
# Replace make's built-in default (FC=f77) with ifx, while still honouring an
# explicit `make FC=...` on the command line or an exported FC in the env.
# (`?=` alone does not override the built-in default, whose origin is 'default'.)
ifeq ($(origin FC),default)
  FC := ifx
endif
FC ?= ifx
# other examples:
# FC = ifort
# FC = ifx
# FC = amdflang
# FC = gfortran

# --- Compiler family ---------------------------------------------------------
# Every flag/library block below keys off FCFAM, never off FC itself: FC is
# often an absolute path (FC=/opt/AMD/aocc-compiler-5.2.0/bin/flang) or a
# wrapper script (FC=.../aoccflang), and matching those against the literal
# string "flang" used to fail silently -- the build then got no -march, no
# -fopenmp and, worst of all, an empty LIBS, i.e. a libfmod.so with unresolved
# BLAS symbols that only blows up at import time. Match on the basename.
# Covered names: ifort/ifx/gfortran, and anything ending in "flang"
# (flang, amdflang, aoccflang, flang-new).
FCNAME := $(notdir $(FC))
ifneq ($(filter %flang %flang-new,$(FCNAME)),)
  FCFAM := flang
else ifeq ($(FCNAME),ifx)
  FCFAM := ifx
else ifeq ($(FCNAME),ifort)
  FCFAM := ifort
else ifeq ($(FCNAME),gfortran)
  FCFAM := gfortran
else
  $(error Unrecognised Fortran compiler FC='$(FC)' (basename '$(FCNAME)'). \
    Use ifx, ifort, gfortran, or a driver whose name ends in "flang". \
    Note FC is also picked up from the environment -- `make FC=ifx` overrides it.)
endif

# C++ compiler for fft_wrapper.cpp: auto-selected to match the FC family.
# (Unused on the MKL FFT path, which needs no C++ wrapper -- see below.)
# Explicit `make CXX=...` on the command line takes priority;
# environment CXX is ignored so FC always picks the matching compiler.
ifeq ($(FCFAM),ifort)
  _CXX_DEFAULT := icpc
else ifeq ($(FCFAM),ifx)
  _CXX_DEFAULT := icpx
else ifeq ($(FCFAM),flang)
  # Prefer the C++ driver sitting next to the Fortran driver, so a path- or
  # wrapper-style FC keeps both halves of the toolchain in the same install.
  _CXX_SIBLING := $(wildcard $(dir $(FC))aoccclang++ $(dir $(FC))clang++)
  _CXX_DEFAULT := $(if $(_CXX_SIBLING),$(firstword $(_CXX_SIBLING)),clang++)
else
  _CXX_DEFAULT := g++
endif

ifneq ($(origin CXX),command line)
  CXX := $(_CXX_DEFAULT)
endif
CXXFLAGS := -O3 -march=$(ARCH) -mtune=$(ARCH) -funroll-loops -fPIC -std=c++17

# Scalar libraries (MKL/AOCL/OpenBLAS); override via SL variable, e.g. SL=OpenBLAS.
# Default is compiler-aware: flang/amdflang (AMD's compiler) pair with AMD's AOCL
# (BLIS/libFLAME, provided via `module load AMD/AOCL`); everything else uses MKL.
ifeq ($(FCFAM),flang)
  SL ?= AOCL
endif
SL ?= MKL
# SL = AOCL
# SL = OpenBLAS

# --- AOCL location -----------------------------------------------------------
# `module load AMD/AOCL` exports AOCL_ROOT. Link against it explicitly rather
# than relying on LIBRARY_PATH, so the build also works when only AOCL_ROOT is
# set, and -rpath makes the resulting .so loadable from Python without the
# module being loaded.
ifneq ($(AOCL_ROOT),)
  AOCL_LIBDIR := $(AOCL_ROOT)/lib
  AOCL_LDFLAGS := -L$(AOCL_LIBDIR) -Wl,-rpath,$(AOCL_LIBDIR)
endif

SRCDIR := src
# Per-compiler, per-arch object dir: objects/modules from different compilers
# are not interchangeable, and znver3/znver4 builds (see `make dispatch`) must
# never mix either.
OBJDIR := obj/$(FCNAME)-$(ARCH)
FLIBSO := libfmod.so
# automatically track all source files in src/
SOURCES := $(wildcard $(SRCDIR)/*.f90)
# create object list from sources; fmod.f90 provides constant module
OBJS := $(patsubst $(SRCDIR)/%.f90,$(OBJDIR)/%.o,$(SOURCES))
# ensure the object for fmod.f90 (constant module) is first
OBJS := $(OBJDIR)/fmod.o $(filter-out $(OBJDIR)/fmod.o,$(OBJS))
# FFT backend selection (the sources call the legacy dfftw_* Fortran entry points):
#   MKL provides its own FFTW3-compatible Fortran wrappers, under the ISSL --
#   build them from $MKLROOT/interfaces/fftw3xf/ first.
#   Everything else uses the bundled pocketfft (BSD-3-Clause) via fft_wrapper.cpp.
# ifort/ifx always use MKL; gfortran uses MKL only when SL=MKL.
#
# NOTE: AOCL does ship a libfftw3.so that already exports dfftw_plan_dft_ &c.,
# so linking -lfftw3 would drop the C++ wrapper entirely -- but that is upstream
# FFTW, which is GPL, and linking it would put libfmod.so as a whole under the
# GPL. Deliberately not used; keep the FFT backend to MKL or pocketfft.
ifeq ($(FCFAM),ifort)
  USE_MKL_FFT := 1
else ifeq ($(FCFAM),ifx)
  USE_MKL_FFT := 1
else ifeq ($(FCFAM),gfortran)
  ifeq ($(SL),MKL)
    USE_MKL_FFT := 1
  endif
endif

ifdef USE_MKL_FFT
  FFT_OBJ  :=
  FFT_DEP  :=
  FFT_LIBS :=
else
  FFT_OBJ  := $(OBJDIR)/fft_wrapper.o
  FFT_DEP  :=
  ifeq ($(FCFAM),flang)
    # ROCm clang++ defaults to GCC 14 headers; point it at the installed GCC.
    # Also pass libstdc++.so by full path since ld.lld skips GCC toolchain dirs.
    GCC_TOOLCHAIN  := $(shell g++ -print-file-name=crtbegin.o | xargs dirname)
    CXXFLAGS       += --gcc-install-dir=$(GCC_TOOLCHAIN)
    FFT_LIBS       := $(shell g++ -print-file-name=libstdc++.so)
  else
    FFT_LIBS := -lstdc++
  endif
endif
# dependency files generated by compiler
DEPFILES := $(OBJS:.o=.d)


FFLAGS_COMMON := -O3 -fPIC
# Keep .mod files in the object directory for every compiler. Both flang
# drivers accept -J (alias of -module-dir) for writing *and* searching, and the
# Intel drivers spell the same thing -module. Leaving them in the working
# directory is not merely untidy: a .mod written by one compiler makes the next
# compiler abort on it ("Reading module ... Unexpected EOF" from gfortran,
# "#7013: not generated by any release of this compiler" from ifx).
# Kept out of FFLAGS_COMMON because that is reused for the link step, where
# flang warns "argument unused during compilation: '-J...'".
ifneq ($(filter $(FCFAM),gfortran flang),)
  MODDIRFLAG := -J$(OBJDIR)
else ifneq ($(filter $(FCFAM),ifort ifx),)
  MODDIRFLAG := -module $(OBJDIR)
endif
FFLAGS_DEBUG := -g -O0

# Helpers for rewriting colon-separated environment path lists (see below).
_EMPTY :=
_SPACE := $(_EMPTY) $(_EMPTY)
_pathjoin = $(subst $(_SPACE),:,$(strip $1))
_pathsplit = $(subst :,$(_SPACE),$1)

# Having several compiler modules loaded at once makes their environment
# variables collide, and the loser is whichever toolchain we are actually
# invoking. Both fixes below are about that, not about the compilers themselves.
ifneq ($(filter $(FCFAM),ifort ifx),)
  # `module load AMD/AOCC` puts /opt/AMD/aocc-compiler-*/include on C_INCLUDE_PATH,
  # which ifx also searches for .mod files. It finds AOCC's omp_lib.mod first and
  # dies with "#7013: This module file was not generated by any release of this
  # compiler", which then cascades into bogus errors about iso_c_binding.
  # The sources contain no INCLUDE statements at all, so dropping foreign
  # toolchain include dirs costs nothing.
  CPATH          := $(call _pathjoin,$(filter-out /opt/AMD/% /opt/rocm%,$(call _pathsplit,$(CPATH))))
  C_INCLUDE_PATH := $(call _pathjoin,$(filter-out /opt/AMD/% /opt/rocm%,$(call _pathsplit,$(C_INCLUDE_PATH))))
  export CPATH C_INCLUDE_PATH
endif

ifeq ($(FCFAM),ifort)
  FFLAGS_COMMON += -mavx2  -qopenmp -heap-arrays
  LIBS := -lmkl_intel_lp64 -lmkl_sequential -lmkl_core $(FFT_LIBS)
endif
ifeq ($(FCFAM),ifx)
  # ifx is LLVM-based: -march=$(ARCH) tunes for the target Zen core and emits
  # vector code with no Intel CPUID dispatch (znver4 -> AVX-512, znver3 -> AVX2).
  # ifx does not accept -mtune=znver*, so -march alone carries the tuning.
  FFLAGS_COMMON += -march=$(ARCH) -funroll-loops \
                   -qopenmp -init:snan -finline-functions -heap-arrays
  LIBS := -lmkl_intel_lp64 -lmkl_sequential -lmkl_core $(FFT_LIBS)
endif
ifeq ($(FCFAM),gfortran)
  # gfortran needs preprocessing enabled when using -MMD
  # also use free form and unlimited line length to avoid truncation errors
  FFLAGS_COMMON += -cpp -ffree-form -ffree-line-length-none \
                   -march=$(ARCH) -mtune=$(ARCH) -funroll-loops \
                   -fopenmp -finline-functions -fimplicit-none \
                   -finit-integer=NaN -finit-real=NaN
  LIBS= -lopenblas $(FFT_LIBS)
  ifeq ($(SL),MKL)
    LIBS := -lmkl_rt $(FFT_LIBS)
  endif
  ifeq ($(SL),AOCL)
    LIBS := -lblis -lflame -lm -lpthread $(FFT_LIBS)
    LDFLAGS_EXTRA := $(AOCL_LDFLAGS)
  endif
endif
# flang / amdflang (ROCm's LLVM flang): inlining is automatic at -O3, and this
# AMD flang has no -finline-functions flag, so it is omitted. BLAS/LAPACK is
# SL-driven; default AOCL (see the compiler-aware SL block above), loaded via
# `module load AMD/AOCL`. The single-threaded libblis (not libblis-mt) is used
# for the same reason MKL uses mkl_sequential: BLAS is called from inside our
# own OpenMP regions, so a threaded BLAS would oversubscribe. MKL likewise uses
# the sequential layer to avoid clashing with flang's own libomp.
ifeq ($(FCFAM),flang)
  FFLAGS_COMMON += -march=$(ARCH) -mtune=$(ARCH) -funroll-loops -fopenmp
  LDFLAGS_EXTRA := $(AOCL_LDFLAGS)
  LIBS := -lblis -lflame -lm -lpthread $(FFT_LIBS)
  ifeq ($(SL),MKL)
    LIBS := -lmkl_intel_lp64 -lmkl_sequential -lmkl_core $(FFT_LIBS)
  endif
  ifeq ($(SL),OpenBLAS)
    LIBS := -lopenblas $(FFT_LIBS)
  endif
  # The AOCC module puts a 32-bit library directory on LIBRARY_PATH, and clang
  # injects every LIBRARY_PATH entry as -L *ahead* of our own -L flags. ld.lld
  # then finds aocc .../lib32/libomp.so first and aborts with
  #   "is incompatible with elf64-x86-64".
  # Drop any lib32 entry -- a 32-bit directory has no role in this link.
  LIBRARY_PATH := $(call _pathjoin,\
                    $(filter-out %/lib32 %/lib32/,$(call _pathsplit,$(LIBRARY_PATH))))
  export LIBRARY_PATH

  # --- AMD LibM ---------------------------------------------------------------
  # ROCm's flang, unlike AOCC's, does not route math intrinsics through AMD's
  # math library on its own. Two independent mechanisms bring it back, and both
  # are needed because they cover different call sites:
  #
  #  * -fveclib=AMDLIBM lets the vectoriser emit amd_vrd4_*/amd_vrd8_* calls for
  #    exp/log/pow/sin/cos/tan/sinh/cosh/erf/... (~5x on an exp-heavy loop --
  #    also faster than AOCC, which does not vectorise exp either). Those names
  #    are unique to AMD LibM, so nothing can interpose them. NOTE there is no
  #    vector tanh in AMD LibM, and get_ffermi is tanh-bound, hence the second
  #    mechanism.
  #  * The scalar names (exp, tanh, ...) are a different matter: libfmod.so is
  #    dlopen'd by ctypes and glibc's libm is already in the process's global
  #    scope, so a plain -lalm always loses the lookup (verified with
  #    LD_DEBUG=bindings). Linking libalm.a *statically* under
  #    -Bsymbolic-functions binds them inside the library instead (~1.2x on
  #    get_ffermi); --exclude-libs keeps them from leaking back out.
  #    libalm.a needs au_cpuid_* from libaoclutils.
  #
  # Cost: results move by <=1 ULP (exp) / <=2 ULP (tanh) against the glibc/MKL
  # build -- physically irrelevant, but no longer bit-identical. `make AMDLIBM=0`
  # opts out when an exact match with the ifx build matters.
  AMDLIBM ?= 1
  _HAS_VECLIB := $(shell $(FC) -fveclib=AMDLIBM --version >/dev/null 2>&1 && echo yes)
  ifeq ($(AMDLIBM)$(_HAS_VECLIB),1yes)
    ifneq ($(wildcard $(AOCL_LIBDIR)/libalm.a),)
      FFLAGS_COMMON += -fveclib=AMDLIBM
      LIBS := -Wl,-Bsymbolic-functions -Wl,--exclude-libs,ALL \
              $(AOCL_LIBDIR)/libalm.a -laoclutils $(LIBS)
    endif
  endif
endif

.PHONY: all clean dispatch

all: $(FLIBSO)

# Runtime CPU dispatch across machines: build BOTH a portable AVX2 library
# (libfmod_znver3.so) and an AVX-512 one (libfmod_znver4.so). _loader.py picks
# libfmod_znver4.so at import when the running CPU has AVX-512 (Zen 4), else the
# znver3 one -- so ONE checkout runs everywhere yet uses AVX-512 only on Zen 4.
# (A single fat binary is not viable here: Intel -ax auto-dispatch gates on a
#  genuine-Intel CPUID and would never take the AVX-512 path on AMD, and Fortran
#  has no per-function multiversioning like C's target_clones.)
dispatch:
	$(MAKE) ARCH=znver3 FLIBSO=libfmod_znver3.so
	$(MAKE) ARCH=znver4 FLIBSO=libfmod_znver4.so
	@echo "Built libfmod_znver3.so (AVX2) and libfmod_znver4.so (AVX-512); _loader.py selects at runtime."

$(OBJDIR):
	mkdir -p $@

# compile each source into object.
# Generate a dependency file only for compilers that accept -MMD; ifx and both
# flang drivers reject it ("unknown argument: '-MMD'"), so they are excluded.
DEPFLAGS = $(if $(filter $(FCFAM),ifx flang),,-MMD -MF $(@:.o=.d))
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS_COMMON) $(MODDIRFLAG) $(if $(DEBUG),$(FFLAGS_DEBUG)) -c $< -o $@ $(DEPFLAGS)

# compile C++ pocketfft wrapper
$(FFT_OBJ): $(SRCDIR)/fft_wrapper.cpp $(SRCDIR)/pocketfft_hdronly.h | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(if $(DEBUG),$(FFLAGS_DEBUG)) -c $< -o $@

# module rule: constant.mod comes from fmod.o
$(OBJDIR)/constant.mod: $(OBJDIR)/fmod.o
	@# module file produced as side effect of compiling fmod.f90

constant.mod: $(OBJDIR)/constant.mod
	@# shortcut to satisfy dependencies that mention unprefixed name

# ensure that every other object waits for fmod.o (order-only)
$(filter-out $(OBJDIR)/fmod.o,$(OBJS)): | $(OBJDIR)/fmod.o

$(FLIBSO): $(OBJS) $(FFT_OBJ) $(FFT_DEP)
	$(FC) $(FFLAGS_COMMON) $(if $(DEBUG),$(FFLAGS_DEBUG)) -shared $(OBJS) $(FFT_OBJ) -o $@ \
		$(LDFLAGS_EXTRA) $(LIBS)

# include generated dependency files (ignore missing initially)
-include $(DEPFILES)

clean:
	rm -rf obj $(FLIBSO) libfmod_znver3.so libfmod_znver4.so *.mod *.d
