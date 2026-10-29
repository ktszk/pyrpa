# ============================================================================
#  pyrpa Fortran kernel -- CMake build
#
#  Mirrors libs/makefile, including its auto-detection:
#
#    cmake -S libs -B build && cmake --build build -j
#        toolchain from what is installed / loaded, CPU targets from the
#        cluster's node inventory (libs/detect_arch.sh) -> one
#        libfmod_<target>.so per distinct node CPU, exactly like `make auto`.
#    cmake -S libs -B build -DARCH=znver4          (Intel: -DARCH=x86-64-v4)
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
option(FFTZ     "use AOCL-FFTZ on the flang path when found (OFF = pocketfft)" ON)
if(DEFINED ENV{FFTZ_ROOT})
  set(_FFTZ_ROOT_DEFAULT "$ENV{FFTZ_ROOT}")
else()
  set(_FFTZ_ROOT_DEFAULT "$ENV{HOME}/opt/aocl-fftz")
endif()
set(FFTZ_ROOT "${_FFTZ_ROOT_DEFAULT}" CACHE PATH "AOCL-FFTZ install prefix")

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

# AOCL-FFTZ (flang only; see the makefile): its FFTW3-compatible libfftw3.so
# provides the dfftw_* entry points, so this path needs no C++ either.  Not used
# with SL=MKL, which brings its own FFTW3 interface.
set(USE_FFTZ FALSE)
if(FCFAM STREQUAL "flang" AND FFTZ AND NOT SL STREQUAL "MKL")
  if(EXISTS "${FFTZ_ROOT}/lib/libaocl_fftz.so")
    set(USE_FFTZ TRUE)
    set(FFTZ_LIBDIR "${FFTZ_ROOT}/lib")
  else()
    message(STATUS "AOCL-FFTZ not found under FFTZ_ROOT=${FFTZ_ROOT}; using pocketfft")
  endif()
endif()

if(NOT CMAKE_Fortran_COMPILER)
  set(CMAKE_Fortran_COMPILER "${FC_COMPILER}")
endif()
if(USE_MKL_FFT OR USE_FFTZ)
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
  # clang++ (AOCC/ROCm) picks the newest GCC it finds for headers and
  # libstdc++, which may be one installed without libstdc++ (e.g. GCC 14
  # alongside a complete GCC 13); project()'s compiler check then dies with
  # "unable to find library -lstdc++".  Point it at the GCC that g++ uses.
  # Must be set before project(), hence CMAKE_CXX_FLAGS_INIT.
  get_filename_component(_CXXNAME "${CMAKE_CXX_COMPILER}" NAME)
  if(_CXXNAME MATCHES "clang")
    execute_process(COMMAND g++ -print-file-name=crtbegin.o
                    OUTPUT_VARIABLE _crtbegin OUTPUT_STRIP_TRAILING_WHITESPACE
                    RESULT_VARIABLE _gxx_rc ERROR_QUIET)
    if(_gxx_rc EQUAL 0 AND IS_ABSOLUTE "${_crtbegin}")
      get_filename_component(_gcc_dir "${_crtbegin}" DIRECTORY)
      string(APPEND CMAKE_CXX_FLAGS_INIT " --gcc-install-dir=${_gcc_dir}")
    endif()
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
    # portable AVX2 level of the build host's vendor, as the makefile's default ARCH
    file(STRINGS /proc/cpuinfo _vendor REGEX "^vendor_id" LIMIT_COUNT 1)
    if(_vendor MATCHES "GenuineIntel")
      set(ARCH_LIST "x86-64-v3")
    else()
      set(ARCH_LIST "znver3")
    endif()
    message(WARNING "detect_arch.sh unavailable; falling back to ARCH_LIST=${ARCH_LIST}")
  endif()
endif()

# -mtune takes CPU names only (an x86-64-vN level tunes for 'generic'); ifort knows
# neither znver* nor the levels: plain ISA switches for those (no Intel CPUID check,
# unlike -x<code>), its own -march for Intel CPU names -- as in the makefile
function(arch_flags arch out_march out_mtune)
  if(arch MATCHES "^x86-64")
    set(${out_mtune} "generic" PARENT_SCOPE)
  else()
    set(${out_mtune} "${arch}" PARENT_SCOPE)
  endif()
  if(FCFAM STREQUAL "ifort")
    if(arch MATCHES "^znver" OR arch STREQUAL "x86-64-v3")
      set(${out_march} "-mavx2" PARENT_SCOPE)
    elseif(arch STREQUAL "x86-64-v4")
      set(${out_march} "-march=skylake-avx512" PARENT_SCOPE)
    elseif(arch STREQUAL "x86-64" OR arch STREQUAL "x86-64-v2")
      set(${out_march} "-msse4.2" PARENT_SCOPE)
    else()
      set(${out_march} "-march=${arch}" PARENT_SCOPE)
    endif()
  else()
    set(${out_march} "-march=${arch}" PARENT_SCOPE)
  endif()
endfunction()

# ---- flags and libraries, per compiler family ------------------------------
set(FFLAGS -O3 -fPIC)
set(CXXFLAGS -O3 -funroll-loops -fPIC -std=c++17)

if(FCFAM STREQUAL "ifort")
  list(APPEND FFLAGS -qopenmp -heap-arrays)   # the ISA switch comes from arch_flags
  # threaded MKL layer (same libiomp5 as -qopenmp): see the note in makefile
  set(BLAS_LIBS -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core)
elseif(FCFAM STREQUAL "ifx")
  # ifx is LLVM-based: -march=<target> emits vector code with no Intel CPUID
  # dispatch, so it runs at full width on AMD as well as Intel.  It rejects -mtune=znver*.
  list(APPEND FFLAGS -funroll-loops -qopenmp -init:snan -finline-functions -heap-arrays)
  # threaded MKL layer (same libiomp5 as -qopenmp): see the note in makefile
  set(BLAS_LIBS -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core)
elseif(FCFAM STREQUAL "gfortran")
  list(APPEND FFLAGS -cpp -ffree-form -ffree-line-length-none -funroll-loops
                     -fopenmp -finline-functions -fimplicit-none
                     -finit-integer=NaN -finit-real=NaN)
  if(SL STREQUAL "MKL")
    # GNU threading layer: gfortran's libgomp, as in the makefile
    set(BLAS_LIBS -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core)
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

# AOCL-FFTZ library lookup: AOCL's lib dir holds AMD-FFTW's libfftw3.so (GPL)
# under the very same soname, and `module load AMD/AOCL` puts it on
# LD_LIBRARY_PATH, which beats a modern DT_RUNPATH.  So the FFTZ dir goes
# *before* AOCL's in -L/-rpath, and --disable-new-dtags emits an old-style
# DT_RPATH, which the loader searches before LD_LIBRARY_PATH.
if(USE_FFTZ)
  list(APPEND LINK_EXTRA "-L${FFTZ_LIBDIR}" "-Wl,--disable-new-dtags"
                         "-Wl,-rpath,${FFTZ_LIBDIR}")
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
# needs no C++ wrapper; flang uses AOCL-FFTZ (BSD-3-Clause) when found, and
# everything else the bundled pocketfft.  AOCL's own libfftw3 is AMD-FFTW,
# i.e. upstream FFTW (GPL), and is deliberately not linked.
if(USE_MKL_FFT)
  set(FFT_SRC "")
  set(FFT_LIBS "")
elseif(USE_FFTZ)
  set(FFT_SRC "")
  # Full path, not -lfftw3: LIBRARY_PATH may put AMD-FFTW (GPL) ahead of FFTZ.
  set(FFT_LIBS "${FFTZ_LIBDIR}/libfftw3.so" -laocl_fftz)
else()
  set(FFT_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/fft_wrapper.cpp")
  if(FCFAM STREQUAL "flang")
    # The C++ side already carries --gcc-install-dir (set before project());
    # the link is driven by flang, so pass libstdc++ by full path since ld.lld
    # skips GCC toolchain dirs.
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
  arch_flags("${arch}" _fmarch _mtune)
  target_compile_options(${outname} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${FFLAGS};${_fmarch}>
    $<$<COMPILE_LANGUAGE:CXX>:${CXXFLAGS};-march=${arch};-mtune=${_mtune}>)
  # ifx has no -mtune=znver*; every other family takes it
  if(NOT FCFAM STREQUAL "ifx" AND NOT FCFAM STREQUAL "ifort")
    target_compile_options(${outname} PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-mtune=${_mtune}>)
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
elseif(USE_FFTZ)
  message(STATUS "FFT backend      : AOCL-FFTZ (${FFTZ_LIBDIR})")
else()
  message(STATUS "FFT backend      : pocketfft (bundled)")
endif()
message(STATUS "Fortran flags    : ${_fflags} -march=<target>")
message(STATUS "Targets          : ${ARCH_LIST}")
message(STATUS "Output           : ${CMAKE_CURRENT_SOURCE_DIR}/lib<name>.so")
message(STATUS "=================================")
