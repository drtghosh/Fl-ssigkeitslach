cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    	test_main.cpp
	sample_tests.cpp
	test_kernel.cpp
	test_marching_cubes.cpp
	test_nSearch.cpp
	test_sph.cpp)
    
add_executable(tests ${SOURCE_FILES})
add_dependencies(tests SPH)
target_link_libraries(tests SPH)
