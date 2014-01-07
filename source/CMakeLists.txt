# Notice name prefix of this variable, set by CMake according
# to value given with "project()" in the root CMakeLists.txt.
include_directories(../include)

add_library(Tmass Tmass.cpp)
add_library(Tshape Tshape.cpp)
add_library(utility utility.cpp)
add_library(singconfig grid_config.cpp)

add_executable(Mesh mesh.cpp)
target_link_libraries(Mesh Tmass Tshape singconfig utility)

#----example code ----
#include_directories(${MyProject_SOURCE_DIR}/MainProject/inc)
# I assume you want to use LibProject as a library in MainProject.
#include_directories(${MyProject_SOURCE_DIR}/LibProject/inc)
#link_directories(${MyProject_SOURCE_DIR}/LibProject/src)
#add_executable(MainProject main.cpp)
#target_link_libraries(MainProject LibProject)