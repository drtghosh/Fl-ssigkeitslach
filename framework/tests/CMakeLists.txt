cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    	test_main.cpp
	test_scene.cpp
	sample_tests.cpp)
    
add_executable(tests ${SOURCE_FILES})
add_dependencies(tests SPH)
target_link_libraries(tests SPH)
add_dependencies(tests PBF)
target_link_libraries(tests PBF)
