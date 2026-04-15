file(GLOB SOURCES "*.f90")

# Ensure fmod.f90 (module file) is first in the list
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/fmod.f90")
list(INSERT SOURCES 0 "${CMAKE_CURRENT_SOURCE_DIR}/fmod.f90")

add_library(src OBJECT ${SOURCES})

# Set up module dependencies
# All other objects depend on fmod.o (order-only dependency)
set_property(SOURCE ${SOURCES} APPEND PROPERTY
    OBJECT_DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/src.dir/fmod.f90.o"
)
