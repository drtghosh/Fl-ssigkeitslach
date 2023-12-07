# External dependencies
## OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

## Eigen
set( EIGEN3_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/extern/Eigen )
if( NOT EIGEN3_INCLUDE_DIR )
    message( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
endif()

set(SOURCE_FILES
	pbf.h
	pbf.cpp
	params.h
	config.h)
    
# Create a library version so it can be linked with other modules
add_library(PBF ${SOURCE_FILES})
add_dependencies(PBF kernel)
target_link_libraries(PBF kernel)
add_dependencies(PBF geometry)
target_link_libraries(PBF geometry)
add_dependencies(PBF mcubes)
target_link_libraries(PBF mcubes)
target_link_libraries(PBF CompactNSearch)
target_include_directories(PBF PUBLIC ${EIGEN3_INCLUDE_DIR} )
target_include_directories(PBF PUBLIC ${PROJECT_SOURCE_DIR}/extern ) # So you can include header files under PBF/extern
