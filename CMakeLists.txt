# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project( MongeAmpereDG)
subdirs(source include)

cmake_minimum_required(VERSION 2.6.2)
# ----------------
 
# Set a default build type for single-configuration
# CMake generators if no build type is set.
IF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE DEBUG)
ENDIF(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)

SET(CMAKE_BUILD_TYPE Debug)

# Select flags.
SET(CMAKE_CXX_FLAGS "-Wall")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g -frounding-math")

message( "We use ${CMAKE_BUILD_TYPE}" )

SITE_NAME(MY_SITE)
IF(MY_SITE STREQUAL "indigo7")
   message( "We are on ${mySite}, use IGPM-specific config." )
   SET(ARPACK_LIBRARY "-L $ENV{HOME}/lib/ARPACK/ -larpack -lgfortran")
ELSE()
   message( "We are on ${mySite}, use general config." )
   SET(ARPACK_LIBRARY "-L /usr/lib/ -larpack -llapack")
ENDIF()

include_directories(${IGPM_T2_LIB_INCLUDES} ${EIGEN3_INCLUDE_DIR} ${UMFPACK_INCLUDES})

find_package(CGAL QUIET COMPONENTS Core )
#find_package(BLAS REQUIRED)
#find_package(IGPM_T2_LIB REQUIRED)
#find_package(Eigen3 REQUIRED)
#find_package(Umfpack REQUIRED)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

    include_directories (BEFORE "include")


else()
  
    message(STATUS "This program requires the CGAL library, and will not be compiled.")
  
endif()


set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
