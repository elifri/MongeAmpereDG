project(dune-mongeampere CXX)


cmake_minimum_required(VERSION 2.8.6)

IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE DEBUG)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET(CMAKE_BUILD_TYPE RELEASE)

# Select flags.
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -fmax-errors=3  -ftemplate-backtrace-limit=0 -fopenmp -Wno-error=strict-overflow")
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
#set(METHOD_OF_ELLIPSOIDS_PATH "$ENV{HOME}/workspace/Method\ of\ Ellipsoids")

#unset(LIBFFO_CORE_INCLUDES CACHE)
#unset(LIBFFO_CORE_LIBRARIES CACHE)
#unset(Boost_FOUND CACHE)
find_package(Umfpack REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(adolc)
find_package(MethodOfEllipsoids)
find_package(CImg)
find_package(HDF5 REQUIRED)
find_package(libffo_core REQUIRED)


include_directories(${LIBFFO_CORE_INCLUDES} ${METHOD_OF_ELLIPSOIDS_INCLUDES} ${CIMG_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${UMFPACK_INCLUDES} ${Boost_INCLUDE_DIRS})
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
add_subdirectory("m4")
add_subdirectory("dune")
add_subdirectory("doc")
add_subdirectory("cmake/modules")

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)

message("compiling mode: ${CMAKE_BUILD_TYPE}")
message("cmake flags mode: ${CMAKE_CXX_FLAGS}")
message("cmake release flags mode: ${CMAKE_CXX_FLAGS_RELEASE}")
