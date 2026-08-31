# Source list for the parent CMakeLists.  fmod.f90 first: it defines the
# `constant` module every other file uses.  (CMake's Fortran module dependency
# scanner handles the rest, so no explicit OBJECT_DEPENDS is needed.)
file(GLOB _SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.f90")
list(REMOVE_ITEM _SRC "${CMAKE_CURRENT_SOURCE_DIR}/fmod.f90")
list(INSERT _SRC 0 "${CMAKE_CURRENT_SOURCE_DIR}/fmod.f90")
set(FMOD_SOURCES ${_SRC} PARENT_SCOPE)
