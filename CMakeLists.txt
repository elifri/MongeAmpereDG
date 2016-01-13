project(dune-mongeampere CXX)


cmake_minimum_required(VERSION 2.8.6)

IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE DEBUG)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET(CMAKE_BUILD_TYPE RELEASE)

# Select flags.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fmax-errors=3  -ftemplate-backtrace-limit=0 -fopenmp")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")



message("compiling mode: ${CMAKE_BUILD_TYPE}")
message("cmake flags mode: ${CMAKE_CXX_FLAGS}")
message("cmake release flags mode: ${CMAKE_CXX_FLAGS_RELEASE}")


SITE_NAME(MY_SITE)
IF(MY_SITE STREQUAL "indigo7")
   message( "We are on ${MY_SITE}, use dogleg method." )
   SET(SNES_LIBRARY Dogleg)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_DOGLEG")
   SET(ADOLC_PATH $ENV{HOME}/lib/adolc)
   unset(Boost_FOUND CACHE)
   SET(Boost_LIBRARIES /usr/local/gcc-5.2.0/lib/libboost_program_options.so)
   SET(BoostPO_LIBRARIES ${Boost_LIBRARIES})
ELSE()
#   message( "We are on ${MY_SITE}, use petsc snes." )
#   SET(SNES_LIBRARY petsc)
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PETSC")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_DOGLEG")
   message( "We are on ${MY_SITE}, use dogleg method." )
   SET(SNES_LIBRARY Dogleg)
   SET(ADOLC_PATH $ENV{HOME}/lib/adolcRapunzel)
   SET(Boost_LIBRARIES boost_program_options)
   SET(BoostPO_LIBRARIES ${Boost_LIBRARIES})
   SET(ARPACK_LIBRARY "-lcholmod -lsuitesparseconfig -lcolamd -lmetis -lccolamd -lcamd")
ENDIF()


if(NOT (dune-common_DIR OR dune-common_ROOT OR
      "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules"
  ${dune-common_MODULE_PATH})

#include the dune macros
include(DuneMacros)

#fieser hack ... 
set(EIGEN3_INCLUDE_DIR "/home/data/friebel/workspace/eigen")
set(IGPM_T2_LIB_INCLUDES "/home/data/friebel/workspace/igpm_t2_lib/include")
set(METHOD_OF_ELLIPSOIDS_PATH "$ENV{HOME}/workspace/Method\ of\ Ellipsoids")
set(CIMG_PATH "$ENV{HOME}/workspace/CImg-1.6.0")

#unset(CIMG_FOUND CACHE)
find_package(Umfpack REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(adolc REQUIRED)
find_package(MethodOfEllipsoids)
find_package(CImg)


include_directories(${METHOD_OF_ELLIPSOIDS_INCLUDES} ${CIMG_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${UMFPACK_INCLUDES} ${IGPM_T2_LIB_INCLUDES} ${ADOLC_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
message("found boost ${Boost_LIBRARIES}")


# start a dune project with information from dune.module
dune_project()
message("found boost ${Boost_LIBRARIES}")

add_subdirectory("src")
add_subdirectory("Testing")
add_subdirectory("m4")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("cmake/modules")


message("found boost ${Boost_LIBRARIES}")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)


message("found boost ${Boost_LIBRARIES}")

message("compiling mode: ${CMAKE_BUILD_TYPE}")
message("cmake flags mode: ${CMAKE_CXX_FLAGS}")
message("cmake release flags mode: ${CMAKE_CXX_FLAGS_RELEASE}")
