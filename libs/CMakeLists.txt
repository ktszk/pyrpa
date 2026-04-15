CMAKE_MINIMUM_REQUIRED (VERSION 3.10)
project(pyrpa Fortran CXX)

# ===== User Configuration Options =====
set(FC_COMPILER "flang" CACHE STRING "Fortran compiler: ifort, ifx, flang, amdflang, gfortran")
set(SL "MKL" CACHE STRING "Scalar Libraries: MKL, AOCL, OpenBLAS")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type: Release, Debug")

# ===== Output Settings =====
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_OBJECT_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/obj")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/obj")

# ===== Fortran and C++ Compiler Selection =====
set(COMPILER_FAMILY "")

if(FC_COMPILER STREQUAL "ifort")
  if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "ifort")
  endif()
  if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "icpc")
  endif()
  set(COMPILER_FAMILY "intel")

elseif(FC_COMPILER STREQUAL "ifx")
  if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "ifx")
  endif()
  if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "icpx")
  endif()
  set(COMPILER_FAMILY "intel_new")

elseif(FC_COMPILER STREQUAL "gfortran")
  if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "gfortran")
  endif()
  if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "g++")
  endif()
  set(COMPILER_FAMILY "gnu")

elseif(FC_COMPILER STREQUAL "flang")
  if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "flang")
  endif()
  if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "clang++")
  endif()
  set(COMPILER_FAMILY "flang")

elseif(FC_COMPILER STREQUAL "amdflang")
  if(NOT CMAKE_Fortran_COMPILER)
    set(CMAKE_Fortran_COMPILER "amdflang")
  endif()
  if(NOT CMAKE_CXX_COMPILER)
    set(CMAKE_CXX_COMPILER "clang++")
  endif()
  set(COMPILER_FAMILY "amdflang")

else()
  message(FATAL_ERROR "Unknown compiler: ${FC_COMPILER}. Valid options: ifort, ifx, gfortran, flang, amdflang")
endif()

# Enable language support
enable_language(Fortran)
enable_language(CXX)

message(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Compiler Family: ${COMPILER_FAMILY}")
message(STATUS "Scalar Library: ${SL}")

# ===== Compiler-specific Flags =====
set(FFLAGS_COMMON "-O2 -fPIC")
set(FFLAGS_DEBUG "-g -O0")
set(CXXFLAGS "-O2 -fPIC -std=c++17")
set(FFT_BACKEND_LIBS "")

# ===== Intel Compilers (ifort, ifx) =====
if(COMPILER_FAMILY STREQUAL "intel" OR COMPILER_FAMILY STREQUAL "intel_new")
  list(APPEND FFLAGS_COMMON "-mavx2" "-qopenmp" "-heap-arrays")
  if(COMPILER_FAMILY STREQUAL "intel_new")
    list(APPEND FFLAGS_COMMON "-init:snan" "-finline-functions")
  endif()
  set(USE_MKL_FFT TRUE)
  set(BLAS_LIBS "-lmkl_intel_lp64" "-lmkl_sequential" "-lmkl_core")

# ===== GNU Compiler (gfortran) =====
elseif(COMPILER_FAMILY STREQUAL "gnu")
  list(APPEND FFLAGS_COMMON "-cpp" "-ffree-form" "-ffree-line-length-none"
                            "-mavx2" "-fopenmp" "-finline-functions" "-fimplicit-none"
                            "-finit-integer=NaN" "-finit-real=NaN")
  list(APPEND FFLAGS_COMMON "-J${CMAKE_Fortran_MODULE_DIRECTORY}")

  # Select BLAS library based on SL variable
  if(SL STREQUAL "MKL")
    set(USE_MKL_FFT TRUE)
    set(BLAS_LIBS "-lmkl_rt")
  elseif(SL STREQUAL "AOCL")
    set(USE_MKL_FFT FALSE)
    set(BLAS_LIBS "-lblis" "-lflame" "-lm" "-lpthread")
  else()  # OpenBLAS (default)
    set(USE_MKL_FFT FALSE)
    set(BLAS_LIBS "-lopenblas")
  endif()

# ===== Flang/AMD Flang (pocketfft backend) =====
elseif(COMPILER_FAMILY STREQUAL "flang" OR COMPILER_FAMILY STREQUAL "amdflang")
  list(APPEND FFLAGS_COMMON "-march=znver3" "-fopenmp" "-finline-functions")
  set(USE_MKL_FFT FALSE)
  set(BLAS_LIBS "-lblis" "-lflame" "-lm" "-lpthread")

  # For AMD clang++, find GCC toolchain
  if(COMPILER_FAMILY STREQUAL "amdflang")
    execute_process(COMMAND g++ -print-file-name=crtbegin.o
                    OUTPUT_VARIABLE GCC_CRTBEGIN OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "/[^/]*$" "" GCC_TOOLCHAIN "${GCC_CRTBEGIN}")
    list(APPEND CXXFLAGS "--gcc-install-dir=${GCC_TOOLCHAIN}")

    execute_process(COMMAND g++ -print-file-name=libstdc++.so
                    OUTPUT_VARIABLE STDCXX_LIB OUTPUT_STRIP_TRAILING_WHITESPACE)
    set(FFT_BACKEND_LIBS "${STDCXX_LIB}")
  else()
    set(FFT_BACKEND_LIBS "-lstdc++")
  endif()

else()
  message(FATAL_ERROR "Compiler family not recognized. Please file a bug report.")
endif()

set(FFT_LIBS "${FFT_BACKEND_LIBS}")

# ===== FFT Backend Configuration =====
if(USE_MKL_FFT)
  if(NOT DEFINED MKLROOT)
    message(FATAL_ERROR "MKLROOT is not set. Load Intel oneAPI MKL module before running CMake.")
  endif()

  set(FFTW3XF_INSTALL "${CMAKE_BINARY_DIR}/obj/fftw3xf")

  # Determine FFTW3 Fortran wrapper compiler
  if(COMPILER_FAMILY STREQUAL "gnu")
    set(FFTW3XF_COMPILER "gnu")
  else()
    set(FFTW3XF_COMPILER "intel")
  endif()

  # Custom target to build MKL FFTW3 Fortran wrapper
  add_custom_target(fftw3xf_build ALL
    COMMAND mkdir -p ${FFTW3XF_INSTALL}
    COMMAND $(MAKE) -C ${MKLROOT}/share/mkl/interfaces/fftw3xf libintel64
            MKLROOT=${MKLROOT}
            INSTALL_DIR=${FFTW3XF_INSTALL}
            compiler=${FFTW3XF_COMPILER}
            INSTALL_LIBNAME=libfftw3xf_intel.a
    COMMENT "Building MKL FFTW3 Fortran wrapper"
  )

  list(APPEND FFT_LIBS "-L${FFTW3XF_INSTALL}" "-lfftw3xf_intel")
  set(FFT_SOURCE "")
  set(FFT_NEEDS_BUILD "yes")

else()
  # Use bundled pocketfft wrapper
  set(FFT_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/fft_wrapper.cpp")
  set(FFT_NEEDS_BUILD "no")
endif()

# ===== Add Subdirectories =====
add_subdirectory(src)

# ===== Build FFT Wrapper Object (if not using MKL) =====
if(FFT_SOURCE)
  add_library(fft_wrapper OBJECT ${FFT_SOURCE})
  target_compile_options(fft_wrapper PRIVATE ${CXXFLAGS})
  set(FFT_OBJ $<TARGET_OBJECTS:fft_wrapper>)
  set(FFT_DEP "")
else()
  set(FFT_OBJ "")
  set(FFT_DEP "fftw3xf_build")
endif()

add_library(fmod SHARED $<TARGET_OBJECTS:src> ${FFT_OBJ})

if(FFT_DEP)
  add_dependencies(fmod ${FFT_DEP})
endif()

# ===== Apply Compiler and Linker Flags =====
string(REPLACE ";" " " FFLAGS_COMMON_STR "${FFLAGS_COMMON}")
string(REPLACE ";" " " FFLAGS_DEBUG_STR "${FFLAGS_DEBUG}")
string(REPLACE ";" " " BLAS_LIBS_STR "${BLAS_LIBS}")
string(REPLACE ";" " " FFT_LIBS_STR "${FFT_LIBS}")

set(CMAKE_Fortran_FLAGS_RELEASE "${FFLAGS_COMMON_STR}")
set(CMAKE_Fortran_FLAGS_DEBUG "${FFLAGS_COMMON_STR} ${FFLAGS_DEBUG_STR}")
set(CMAKE_Fortran_FLAGS "${FFLAGS_COMMON_STR}" CACHE STRING "" FORCE)

set(CMAKE_CXX_FLAGS_RELEASE "${CXXFLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CXXFLAGS} -g")

target_link_libraries(fmod PRIVATE ${BLAS_LIBS_STR} ${FFT_LIBS_STR})

# ===== Setup for Python ctypes =====
set_target_properties(fmod PROPERTIES
  PREFIX ""
  SUFFIX ".so"
)

# ===== Print Build Summary =====
message(STATUS "=== Build Configuration ===")
message(STATUS "Compiler Family: ${COMPILER_FAMILY}")
message(STATUS "BLAS Libraries: ${BLAS_LIBS_STR}")
if(USE_MKL_FFT)
  message(STATUS "FFT Backend: MKL FFTW3 (MKLROOT=${MKLROOT})")
else()
  message(STATUS "FFT Backend: pocketfft (bundled)")
endif()
message(STATUS "Fortran Flags: ${FFLAGS_COMMON_STR}")
message(STATUS "Output Library: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/fmod.so")
message(STATUS "==========================")
