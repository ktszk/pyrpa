CMAKE_MINIMUM_REQUIRED (VERSION 3.5)
project (pyrpa Fortran)
add_subdirectory(src)
add_library(fmod SHARED $<TARGET_OBJECTS:src>)

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  target_compile_options(fmod PRIVATE 
      $<$<CONFIG:Release>:-O2 -mavx2 -qopenmp -IPO -fPIC -heap-arrays>
      $<$<CONFIG:Debug>:-g -traceback>)
  target_link_libraries(fmod PRIVATE mkl_sequential mkl_core iomp5)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(fmod PRIVATE 
      $<$<CONFIG:Release>:-O2 -mavx2 -fopenmp -finline-functions -fimplicit-none -finit-integer=NaN -finit-real=NaN -fPIC>
      $<$<CONFIG:Debug>:-g>)
  target_link_libraries(fmod PRIVATE mkl_sequential mkl_core iomp5)
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "flang")
  target_compile_options(fmod PRIVATE 
      $<$<CONFIG:Release>:-O2 -march=znver3 -fopenmp -finline-functions -fPIC>)
  target_link_libraries(fmod PRIVATE blis flame fftw3 m pthread)
endif()

if(MKL_FOUND)
  target_link_libraries(fmod PRIVATE MKL::MKL)
endif()
