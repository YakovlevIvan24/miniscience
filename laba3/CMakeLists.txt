cmake_minimum_required(VERSION 3.5)

set(PROJECT_NAME demo)
project(${PROJECT_NAME})
list (APPEND CMAKE_MODULE_PATH "${MyProject_SOURCE_DIR}/cmake-modules")
# Get DOLFIN configuration data (DOLFINConfig.cmake must be in DOLFIN_CMAKE_CONFIG_PATH)
# set(DOLFIN_CMAKE_CONFIG_PATH "/home/user/me/dolfin/build/dolfin/")
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

add_executable(poisson poisson.cpp poisson.cpp)
target_link_libraries(poisson dolfin)

