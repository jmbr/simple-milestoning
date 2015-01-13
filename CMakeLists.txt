cmake_minimum_required(VERSION 2.8)

project(simple-milestoning CXX)
enable_language(CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CXX_WARN_FLAGS "-W -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-g -std=c++11 ${CXX_WARN_FLAGS}")

add_executable(milestoning milestoning.cpp integrator.cpp milestones.cpp)
target_link_libraries(milestoning ${ARMADILLO_LIBRARIES})
