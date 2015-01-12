cmake_minimum_required(VERSION 2.8)

project(torus CXX)
enable_language(CXX)

find_package(Armadillo REQUIRED PATHS /workspace/jmb)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set(CXX_WARN_FLAGS "-W -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-g -std=c++11 ${CXX_WARN_FLAGS} ${OpenMP_CXX_FLAGS}")

add_executable(milestoning milestoning.cpp integrator.cpp milestones.cpp)
target_link_libraries(milestoning ${ARMADILLO_LIBRARIES})
