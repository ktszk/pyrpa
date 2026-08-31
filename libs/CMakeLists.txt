# ============================================================================
#  pyrpa Fortran kernel -- CMake build
#
#  Mirrors libs/makefile, including its auto-detection:
#
#    cmake -S libs -B build && cmake --build build -j
#        toolchain from what is installed / loaded, CPU targets from the
#        cluster's node inventory (libs/detect_arch.sh) -> one
#        libfmod_<target>.so per distinct node CPU, exactly like `make auto`.
#    cmake -S libs -B build -DARCH=znver4
#        a single libfmod.so for one target, like `make ARCH=znver4`.
#    cmake -S libs -B build -DFC_COMPILER=gfortran -DSL=MKL -DARCH_LIST="znver3;znver4"
#        everything pinned by hand.
#
#  Libraries are written into libs/ (not the build tree) because
#  libs/flibs/_loader.py loads them from there, picking at import the highest
#  ISA the running node can execute.
#
#  NOTE the makefile remains the reference build.  The exotic AMD paths it
#  carries (CPATH sanitising when several compiler modules are loaded, dropping
#  a lib32 entry from LIBRARY_PATH under AOCC) are environment fix-ups that make
#  applies to its own child processes; with flang+AOCC modules prefer `make`.
# ============================================================================
cmake_minimum_required(VERSION 3.18)

# ---- user knobs (empty = auto) ---------------------------------------------
set(FC_COMPILER "" CACHE STRING "Fortran compiler: ifx, ifort, flang, amdflang, gfortran (empty = first found)")
set(SL          "" CACHE STRING "Scalar library: MKL, AOCL, OpenBLAS (empty = from MKLROOT/AOCL_ROOT)")
set(ARCH        "" CACHE STRING "single -march target; builds libfmod.so (empty = use ARCH_LIST)")
set(ARCH_LIST   "" CACHE STRING "-march targets to build (empty = detect_arch.sh)")
option(AMDLIBM  "link AMD LibM on the flang path" ON)

# ---- compiler auto-selection (MUST precede project(): enable_language pins it)
if(NOT FC_COMPILER)
  foreach(_c ifx ifort flang amdflang gfortran)
    find_program(_FC_FOUND ${_c})
    if(_FC_FOUND)
      set(FC_COMPILER ${_c})
      break()
    endif()
    unset(_FC_FOUND CACHE)
  endforeach()
  unset(_FC_FOUND CACHE)
endif()
if(NOT FC_COMPILER)
  message(FATAL_ERROR "No Fortran compiler found (looked for ifx, ifort, flang, amdflang, gfortran)")
endif()

# Family, matched on the basename like the makefile: FC may be an absolute path
# or a wrapper whose name merely ends in "flang".
get_filename_component(_FCNAME "${FC_COMPILER}" NAME)
if(_FCNAME MATCHES "flang")
  set(FCFAM flang)
  set(_CXX_DEFAULT clang++)
elseif(_FCNAME STREQUAL "ifx")
  set(FCFAM ifx)
  set(_CXX_DEFAULT icpx)
elseif(_FCNAME STREQUAL "ifort")
  set(FCFAM ifort)
  set(_CXX_DEFAULT icpc)
elseif(_FCNAME STREQUAL "gfortran")
  set(FCFAM gfortran)
  set(_CXX_DEFAULT g++)
else()
  message(FATAL_ERROR "Unrecognised Fortran compiler '${FC_COMPILER}'. Use ifx, ifort, "
                      "gfortran, or a driver whose name ends in \"flang\".")
endif()

# ---- scalar library: follow the modules actually loaded ---------------------
# Decided before project() because it also decides the FFT backend, and hence
# whether a C++ compiler is needed at all.
if(NOT SL)
  if(FCFAM STREQUAL "flang" AND DEFINED ENV{AOCL_ROOT})
    set(SL AOCL)
  elseif(DEFINED ENV{MKLROOT})
    set(SL MKL)
  elseif(DEFINED ENV{AOCL_ROOT})
    set(SL AOCL)
  else()
    set(SL OpenBLAS)
  endif()
endif()

# MKL ships FFTW3-compatible Fortran entry points (dfftw_*), so that path needs
# no C++ at all -- which matters: an oneAPI install can carry ifx without icpx.
if(FCFAM STREQUAL "ifx" OR FCFAM STREQUAL "ifort"
   OR (FCFAM STREQUAL "gfortran" AND SL STREQUAL "MKL"))
  set(USE_MKL_FFT TRUE)
else()
  set(USE_MKL_FFT FALSE)
endif()

if(NOT CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER "${FC_COMPILER}")
endif()
if(USE_MKL_FFT)
  project(pyrpa Fortran)
else()
  if(NOT CMAKE_CXX_COMPILER)
    # the family's own C++ driver if it is installed, otherwise g++
    find_program(_CXX_FOUND ${_CXX_DEFAULT})
    if(_CXX_FOUND)
      set(CMAKE_CXX_COMPILER "${_CXX_FOUND}")
    else()
      set(CMAKE_CXX_COMPILER "g++")
      message(STATUS "${_CXX_DEFAULT} not found; using g++ for fft_wrapper.cpp")
    endif()
    unset(_CXX_FOUND CACHE)
  endif()
  project(pyrpa Fortran CXX)
endif()

# ---- target CPU list -------------------------------------------------------
# Never -march=native: that targets the BUILD host, and an AVX-512 build started
# on an older node dies with SIGILL.  detect_arch.sh asks the compute nodes.
set(_SINGLE_LIB FALSE)
if(ARCH)
  set(ARCH_LIST "${ARCH}")
  set(_SINGLE_LIB TRUE)          # one target, plain libfmod.so, like `make ARCH=...`
elseif(NOT ARCH_LIST)
  execute_process(COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/detect_arch.sh"
                  OUTPUT_VARIABLE _detected RESULT_VARIABLE _rc
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(_rc EQUAL 0 AND _detected)
    string(REPLACE "\n" ";" ARCH_LIST "${_detected}")
  else()
    set(ARCH_LIST "znver3")      # portable AVX2 baseline
    message(WARNING "detect_arch.sh unavailable; falling back to ARCH_LIST=znver3")
  endif()
endif()

# ---- flags and libraries, per compiler family ------------------------------
set(FFLAGS -O3 -fPIC)
set(CXXFLAGS -O3 -funroll-loops -fPIC -std=c++17)

if(FCFAM STREQUAL "ifort")
  list(APPEND FFLAGS -mavx2 -qopenmp -heap-arrays)
  set(BLAS_LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core)
elseif(FCFAM STREQUAL "ifx")
  # ifx is LLVM-based: -march=<target> emits vector code with no Intel CPUID
  # dispatch, so it runs at full width on AMD.  It rejects -mtune=znver*.
  list(APPEND FFLAGS -funroll-loops -qopenmp -init:snan -finline-functions -heap-arrays)
  set(BLAS_LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core)
elseif(FCFAM STREQUAL "gfortran")
  list(APPEND FFLAGS -cpp -ffree-form -ffree-line-length-none -funroll-loops
                     -fopenmp -finline-functions -fimplicit-none
                     -finit-integer=NaN -finit-real=NaN)
  if(SL STREQUAL "MKL")
    set(BLAS_LIBS -lmkl_rt)
  elseif(SL STREQUAL "AOCL")
    set(BLAS_LIBS -lblis -lflame -lm -lpthread)
  else()
    set(BLAS_LIBS -lopenblas)
  endif()
else()  # flang / amdflang
  list(APPEND FFLAGS -funroll-loops -fopenmp)
  if(SL STREQUAL "MKL")
    set(BLAS_LIBS -lmkl_intel_lp64 -lmkl_sequential -lmkl_core)
  elseif(SL STREQUAL "OpenBLAS")
    set(BLAS_LIBS -lopenblas)
  else()
    # single-threaded libblis (not libblis-mt): BLAS is called from inside our own
    # OpenMP regions, so a threaded BLAS would oversubscribe.
    set(BLAS_LIBS -lblis -lflame -lm -lpthread)
  endif()
endif()

# AOCL location: `module load AMD/AOCL` exports AOCL_ROOT; -rpath keeps the .so
# loadable from Python without the module loaded.
if(DEFINED ENV{AOCL_ROOT})
  set(AOCL_LIBDIR "$ENV{AOCL_ROOT}/lib")
  list(APPEND LINK_EXTRA "-L${AOCL_LIBDIR}" "-Wl,-rpath,${AOCL_LIBDIR}")
endif()

# AMD LibM on the flang path (see the makefile for the two mechanisms and the
# <=2 ULP cost); -DAMDLIBM=OFF opts out.
if(FCFAM STREQUAL "flang" AND AMDLIBM AND DEFINED ENV{AOCL_ROOT}
   AND EXISTS "${AOCL_LIBDIR}/libalm.a")
  execute_process(COMMAND ${CMAKE_Fortran_COMPILER} -fveclib=AMDLIBM --version
                  RESULT_VARIABLE _veclib_rc OUTPUT_QUIET ERROR_QUIET)
  if(_veclib_rc EQUAL 0)
    list(APPEND FFLAGS -fveclib=AMDLIBM)
    list(PREPEND BLAS_LIBS -Wl,-Bsymbolic-functions -Wl,--exclude-libs,ALL
                           "${AOCL_LIBDIR}/libalm.a" -laoclutils)
  endif()
endif()

# ---- FFT backend -----------------------------------------------------------
# MKL ships FFTW3-compatible Fortran entry points (dfftw_*), so the MKL path
# needs no C++ wrapper; everything else uses the bundled pocketfft.  AOCL's
# libfftw3 is upstream FFTW (GPL) and is deliberately not linked.
if(USE_MKL_FFT)
  set(FFT_SRC "")
  set(FFT_LIBS "")
else()
  set(FFT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/fft_wrapper.cpp")
  if(FCFAM STREQUAL "flang")
    # ROCm clang++ defaults to newer GCC headers; point it at the installed GCC,
    # and pass libstdc++ by full path since ld.lld skips GCC toolchain dirs.
    execute_process(COMMAND g++ -print-file-name=crtbegin.o
                    OUTPUT_VARIABLE _crtbegin OUTPUT_STRIP_TRAILING_WHITESPACE)
    get_filename_component(_gcc_dir "${_crtbegin}" DIRECTORY)
    list(APPEND CXXFLAGS "--gcc-install-dir=${_gcc_dir}")
    execute_process(COMMAND g++ -print-file-name=libstdc++.so
                    OUTPUT_VARIABLE FFT_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
  else()
    set(FFT_LIBS -lstdc++)
  endif()
endif()

add_subdirectory(src)          # sets FMOD_SOURCES

# ---- one shared library per target ISA -------------------------------------
function(add_fmod_variant arch outname)
  add_library(${outname} SHARED ${FMOD_SOURCES} ${FFT_SRC})
  set_target_properties(${outname} PROPERTIES
    OUTPUT_NAME "${outname}" PREFIX "lib" SUFFIX ".so"
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod/${arch}"
    # land next to the makefile's output: _loader.py loads from libs/
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  target_compile_options(${outname} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${FFLAGS};-march=${arch}>
    $<$<COMPILE_LANGUAGE:CXX>:${CXXFLAGS};-march=${arch};-mtune=${arch}>)
  # ifx has no -mtune=znver*; every other family takes it
  if(NOT FCFAM STREQUAL "ifx" AND NOT FCFAM STREQUAL "ifort")
    target_compile_options(${outname} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-mtune=${arch}>)
  endif()
  target_link_options(${outname} PRIVATE ${FFLAGS} ${LINK_EXTRA})
  target_link_libraries(${outname} PRIVATE ${BLAS_LIBS} ${FFT_LIBS})
endfunction()

foreach(_arch IN LISTS ARCH_LIST)
  if(_SINGLE_LIB)
    add_fmod_variant(${_arch} fmod)
  else()
    add_fmod_variant(${_arch} fmod_${_arch})
  endif()
endforeach()

# ---- summary ---------------------------------------------------------------
string(REPLACE ";" " " _blas "${BLAS_LIBS}")
string(REPLACE ";" " " _fflags "${FFLAGS}")
message(STATUS "=== pyrpa build configuration ===")
message(STATUS "Fortran compiler : ${CMAKE_Fortran_COMPILER} (family ${FCFAM})")
message(STATUS "Scalar library   : ${SL} -> ${_blas}")
if(USE_MKL_FFT)
  message(STATUS "FFT backend      : MKL (FFTW3 Fortran entry points)")
else()
  message(STATUS "FFT backend      : pocketfft (bundled)")
endif()
message(STATUS "Fortran flags    : ${_fflags} -march=<target>")
message(STATUS "Targets          : ${ARCH_LIST}")
message(STATUS "Output           : ${CMAKE_CURRENT_SOURCE_DIR}/lib<name>.so")
message(STATUS "=================================")
