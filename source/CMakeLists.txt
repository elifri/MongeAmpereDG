# Notice name prefix of this variable, set by CMake according
# to value given with "project()" in the root CMakeLists.txt.
include_directories(../include)
subdirs(test grid)

# util libraries
add_library(utility utility.cpp)
add_library(Eigen_utility Eigen_utility.cpp)
add_library(Callback_utility Callback_utility.cpp)
set(UTILITIES Eigen_utility singconfig utility)

add_library(Plotter Plotter.cpp)
add_library(Convexifier Convexifier.cpp)
add_library(ConvexifyNLP ConvexifyNLP.cpp)

add_library(Tsolver Tsolver_general.cpp Tsolver_special_POISSON.cpp Tsolver_special_AM.cpp)
set(dependecies_Tsolver Convexifier testlibs ConvexifyNLP boundary_handler C0_converter Plotter Tmass Tshape Callback_utility ${UTILITIES} ${IPOPT_LIBRARIES})

add_executable(Mesh mesh.cpp)
target_link_libraries(Mesh Tsolver ${dependecies_Tsolver} ${BLAS_LIBRARIES} ${ARPACK_LIBRARY})

add_executable(Plot plot.cpp)
target_link_libraries(Plot Tsolver ${dependecies_Tsolver})

add_executable(Test unittest.cpp)
target_link_libraries(Test testlibs Tsolver ${dependecies_Tsolver})  

#----example code ----
#include_directories(${MyProject_SOURCE_DIR}/MainProject/inc)
# I assume you want to use LibProject as a library in MainProject.
#include_directories(${MyProject_SOURCE_DIR}/LibProject/inc)
#link_directories(${MyProject_SOURCE_DIR}/LibProject/src)
#add_executable(MainProject main.cpp)
#target_link_libraries(MainProject LibProject)
