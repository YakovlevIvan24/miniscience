cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS "-std=c++11 -g")

project(c++)


add_executable(laba2.1 laba2.1.cpp)
add_executable(laba1.1 laba1.1.cpp)

# Найдём VTK
find_package(VTK REQUIRED)
# ... и подключим её со всеми её компонентами её нативными средствами
include(${VTK_USE_FILE})

# gmsh укажем явно по кусочкам (потому что он был просто распакован руками)
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
    IMPORTED_LOCATION "/home/user/me/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so"
    INTERFACE_INCLUDE_DIRECTORIES "/home/user/me/gmsh-4.11.1-Linux64-sdk/include"
)
set(GMSH_LIBRARY gmsh)

# если бы gmsh был поставлен прилично, то можно было бы примерно вот так
#find_package(gmsh REQUIRED)
#include_directories(${GMSH_INCLUDE_DIR})
#target_link_libraries(tetr3d ${GMSH_LIBRARY})




target_link_libraries(laba2.1 ${VTK_LIBRARIES})
target_link_libraries(laba2.1 ${GMSH_LIBRARY})


target_link_libraries(laba1.1 ${GMSH_LIBRARY})
