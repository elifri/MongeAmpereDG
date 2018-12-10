project(dune-mongeampere CXX)


cmake_minimum_required(VERSION 2.8.6)

IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE DEBUG)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET(CMAKE_BUILD_TYPE DEBUG)

#  Selectflags.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fmax-errors=3  -ftemplate-backtrace-limit=0 -fopenmp -Wno-error=strict-overflow -Wno-deprecated")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")

set(CMAKE_EXE_LINK_DYNAMIC_C_FLAGS "")
set(CMAKE_EXE_LINK_DYNAMIC_CXX_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
set(LINK_SEARCH_START_STATIC "")
set(LINK_SEARCH_END_STATIC "")


message("compiling mode: ${CMAKE_BUILD_TYPE}")
message("cmake flags mode: ${CMAKE_CXX_FLAGS}")
message("cmake release flags mode: ${CMAKE_CXX_FLAGS_RELEASE}")


SITE_NAME(MY_SITE)
IF(MY_SITE STREQUAL "indigo7")
   message( "We are on ${MY_SITE}, use dogleg method." )
   SET(SNES_LIBRARY Dogleg)
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_DOGLEG")
   SET(ADOLC_PATH $ENV{HOME}/lib/adolc)
   #unset(Boost_FOUND CACHE)
#   SET(Boost_LIBRARIES /usr/lib64/libboost_program_options.so)
   SET(Boost_LIBRARIES -lboost_program_options)
   SET(BoostPO_LIBRARIES ${Boost_LIBRARIES})
ELSE()
#   message( "We are on ${MY_SITE}, use petsc snes." )
#   SET(SNES_LIBRARY petsc)
#   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_PETSC")
   SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_DOGLEG")
   message( "We are on ${MY_SITE}, use dogleg method." )
   SET(SNES_LIBRARY Dogleg)
   unset(Boost_FOUND CACHE)
   SET(Boost_LIBRARIES /usr/local/gcc-7.2.0/lib/libboost_program_options.so)
   SET(ADOLC_PATH /home/terra/friebel/workspace/ADOL-C-2.6.3/build)
   #SET(Boost_LIBRARIES boost_program_options)
   SET(BoostPO_LIBRARIES ${Boost_LIBRARIES})
   SET(ARPACK_LIBRARY "-lcholmod -lsuitesparseconfig -lcolamd -lmetis -lccolamd -lcamd")
   SET(LIBNURBS_DIR /home/terra/friebel/workspace/opennurbs_20130711)
   SET(CIMG_PATH /home/terra/friebel/workspace/CImg-1.6.0)
   SET(IPOPT_HOME /home/terra/friebel/workspace/ipopt/Ipopt-3.12.10/build)
   SET(METHOD_OF_ELLIPSOIDS_PATH /home/terra/friebel/workspace/Method\ of\ Ellipsoids)
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
#set(METHOD_OF_ELLIPSOIDS_PATH "$ENV{HOME}/workspace/Method\ of\ Ellipsoids")

#unset(LIBFFO_CORE_INCLUDES CACHE)
#unset(LIBFFO_CORE_LIBRARIES CACHE)
message("eigen3 includes ${EIGEN3_INCLUDE_DIR}")
find_package(Umfpack REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(adolc)
find_package(MethodOfEllipsoids)
find_package(CImg)
find_package(HDF5)
find_package(libffo_core)
find_package(libnurbs REQUIRED)
find_package(IPOPT)

message("include: ${IPOPT_INCLUDE_DIR} and using Boost ${BoostPO_LIBRARIES}")
message("libnurbs includes ${LIBNURBS_INCLUDES}")

include_directories(SYSTEM ${LIBFFO_CORE_INCLUDES} ${LIBNURBS_INCLUDES} ${IPOPT_INCLUDE_DIR} ${METHOD_OF_ELLIPSOIDS_INCLUDES} ${CIMG_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${UMFPACK_INCLUDES} ${Boost_INCLUDE_DIRS})

if (METHOD_OF_ELLIPSOIDS_INCLUDES)
include_directories(${METHOD_OF_ELLIPSOIDS_INCLUDES})
endif (METHOD_OF_ELLIPSOIDS_INCLUDES)
if (ADOLC_FOUND)
include_directories(${ADOLC_INCLUDE_DIRS})
endif (ADOLC_FOUND)
include_directories("include")


# start a dune project with information from dune.module
dune_project()

add_subdirectory("include")
add_subdirectory("src")
add_subdirectory("Testing")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("cmake/modules")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

message("compiling mode: ${CMAKE_BUILD_TYPE}")
message("cmake flags mode: ${CMAKE_CXX_FLAGS}")
message("cmake release flags mode: ${CMAKE_CXX_FLAGS_RELEASE}")
