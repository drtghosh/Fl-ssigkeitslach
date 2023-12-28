cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    	test_main.cpp
	test_emitter
	sample_tests.cpp
	test_kernel.cpp
	test_marching_cubes
	test_marching_cubes_pbf.cpp
	test_marching_cubes_sph.cpp
	test_nSearch.cpp
	test_pbf.cpp
	test_sph.cpp)
    
add_executable(tests ${SOURCE_FILES})
add_dependencies(tests SPH)
target_link_libraries(tests SPH)
add_dependencies(tests PBF)
target_link_libraries(tests PBF)
