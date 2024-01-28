cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
	test_compare.cpp
    	test_main.cpp
	test_emitter_ob.cpp
	sample_tests.cpp
	test_kernel.cpp
	test_nSearch.cpp)
    
add_executable(tests ${SOURCE_FILES})
add_dependencies(tests SPH)
target_link_libraries(tests SPH)
add_dependencies(tests PBF)
target_link_libraries(tests PBF)
