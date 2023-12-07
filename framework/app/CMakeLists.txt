cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)

# Local
set(SOURCE_FILES
    main.cpp)
    
add_executable(app ${SOURCE_FILES})
add_dependencies(app SPH)
target_link_libraries(app SPH)
