# Notice name prefix of this variable, set by CMake according
# to value given with "project()" in the root CMakeLists.txt.
include_directories(../include)

add_library(testlibs init_stiffnessmatrix.cpp test_utils.cpp)