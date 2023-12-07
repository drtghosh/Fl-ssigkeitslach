# Install script for directory: C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files/Fluid")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/extern/CompactNSearch/cmake_install.cmake")
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/app/cmake_install.cmake")
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/geometry/cmake_install.cmake")
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/kernel/cmake_install.cmake")
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/mcubes/cmake_install.cmake")
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/PBF/cmake_install.cmake")
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/SPH/cmake_install.cmake")
  include("C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/tests/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "C:/Users/Acer/Documents/All_of_RWTH/Semester3/Lab/Flüssigkeitslach/framework/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
