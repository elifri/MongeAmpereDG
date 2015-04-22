cmake_minimum_required(VERSION 2.8.6)

IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE DEBUG)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET(CMAKE_BUILD_TYPE DEBUG)


# Select flags.
SET(CMAKE_CXX_FLAGS "-std=c++11 -Wall -fmax-errors=3 -w")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")


SITE_NAME(MY_SITE)
IF(MY_SITE STREQUAL "indigo7")
   message( "We are on ${MY_SITE}, use dogleg method." )
   SET(SNES_LIBRARY Dogleg)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_DOGLEG")
ELSE()
   message( "We are on ${MY_SITE}, use petsc snes." )
   SET(SNES_LIBRARY petsc)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PETSC")
ENDIF()


project(dune-mongeampere CXX)


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
message("CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} \n")

find_package(Umfpack REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(adolc REQUIRED)
#find_package(IGPM_T2_LIB REQUIRED HINTS ${IGPM_T2_LIB_INCLUDE_DIR})

include_directories(${EIGEN3_INCLUDE_DIR} ${UMFPACK_INCLUDES} ${IGPM_T2_LIB_INCLUDES} ${ADOLC_INCLUDES})




# start a dune project with information from dune.module
dune_project()

add_subdirectory("src")
add_subdirectory("m4")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("cmake/modules")
include_directories(${dune-localfunctions_INCLUDE_DIR})

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
