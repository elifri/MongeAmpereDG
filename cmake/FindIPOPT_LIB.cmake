# - Try to find IPOPT
# Once done this will define
#  IPOPT_FOUND - System has IpOpt
#  IPOPT_INCLUDE_DIRS - The IpOpt include directories
#  IPOPT_LIBRARY_DIRS - The library directories needed to use IpOpt
#  IPOPT_LIBRARIES    - The libraries needed to use IpOpt

SET(IPOPT_HOME ${HOME}/lib/ipopt/build)

if (IPOPT_INCLUDE_DIR)
  # in cache already
  SET(IPOPT_FIND_QUIETLY TRUE)
endif (IPOPT_INCLUDE_DIR)

if (WIN32)
   find_path(IPOPT_INCLUDE_DIR NAMES IpNLP.hpp
     PATHS
     "C:\\libs\\Ipopt-3.8.2\\include\\coin"
     ${IPOPT_DIR}/include
   )

   IF(IPOPT_INCLUDE_DIR)
      find_library( IPOPT_LIBRARY_RELEASE 
                    Ipopt
                    PATHS "C:\\libs\\Ipopt-3.8.2\\lib\\win32\\release" )
      find_library( IPOPT_LIBRARY_DEBUG
                    Ipopt
                    PATHS "C:\\libs\\Ipopt-3.8.2\\lib\\win32\\debug" )

      set ( IPOPT_LIBRARY "optimized;${IPOPT_LIBRARY_RELEASE};debug;${IPOPT_LIBRARY_DEBUG}" CACHE  STRING "IPOPT Libraries" )

      SET(IPOPT_FOUND TRUE)
      SET(IPOPT_INCLUDE_DIR ${IPOPT_INCLUDE_DIR})
	  # Todo, set right version depending on build type (debug/release)
	  #GET_FILENAME_COMPONENT( IPOPT_LIBRARY_DIR ${GLEW_LIBRARY} PATH )
    ELSE(IPOPT_INCLUDE_DIR)
      SET(IPOPT_FOUND FALSE)
      SET(IPOPT_INCLUDE_DIR ${IPOPT_INCLUDE_DIR})
    ENDIF(IPOPT_INCLUDE_DIR)

ELSE( WIN32 )
   find_path(IPOPT_INCLUDE_DIR NAMES IpNLP.hpp
     PATHS  "${IPOPT_HOME}/include/coin"
            "/usr/include/coin"
    
   )

   find_library( IPOPT_LIBRARY 
                 ipopt
                 PATHS "${IPOPT_HOME}/lib"
                       "/usr/lib" )   
    
    #wrong config under Debian workaround
    add_definitions( -DHAVE_CSTDDEF )

   
   # set optional path to HSL Solver for dynamic usage
   find_path(IPOPT_HSL_LIBRARY_DIR 
             NAMES libhsl.so
                   libhsl.dylib
             PATHS "$ENV{IPOPT_HSL_LIBRARY_PATH}"
                   "$ENV{HOME}/opt/HSL/lib"
   )

   # find HSL library for fixed linking of solvers   
   find_library( IPOPT_HSL_LIBRARY 
                 coinhsl
                 PATHS "${IPOPT_HOME}/lib"
                       "/usr/lib" )   
   
   
   IF( IPOPT_HSL_LIBRARY_DIR)
     IF( NOT IPOPT_FIND_QUIETLY )
        message ( "IPOPT_HSL_LIBRARY_DIR found at ${IPOPT_HSL_LIBRARY_DIR} ")
     ENDIF()
     set(IPOPT_LIBRARY_DIR ${IPOPT_HSL_LIBRARY_DIR})
     LIST( APPEND IPOPT_LIBRARY_DIRS "${IPOPT_HSL_LIBRARY_DIR}")
   ENDIF(IPOPT_HSL_LIBRARY_DIR)
   
   
   set(IPOPT_INCLUDE_DIRS "${IPOPT_INCLUDE_DIR}" )
   set(IPOPT_LIBRARIES "${IPOPT_LIBRARY}" )
   
   IF(IPOPT_HSL_LIBRARY)
     LIST( APPEND IPOPT_LIBRARIES "${IPOPT_HSL_LIBRARY}")   
   ENDIF(IPOPT_HSL_LIBRARY)

   FIND_LIBRARY(COIN_HSL_LIBRARY coinhsl HINTS  ${IPOPT_HOME}/coin/lib/ ${IPOPT_HOME}/coin/lib/ThirdParty/ ${IPOPT_HOME}/lib/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ /usr/local/lib/coin /usr/local/lib/coin/ThirdParty)	
   FIND_LIBRARY(COIN_MUMPS_LIBRARY  coinmumps HINTS  ${IPOPT_HOME}/coin/lib/ ${IPOPT_HOME}/coin/lib/ThirdParty/ ${IPOPT_HOME}/lib/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ /usr/local/lib/coin /usr/local/lib/coin/ThirdParty)
	
   FIND_LIBRARY(COIN_METIS_LIBRARY coinmetis HINTS  ${IPOPT_HOME}/coin/lib/ ${IPOPT_HOME}/coin/lib/ThirdParty/ ${IPOPT_HOME}/lib/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ /usr/local/lib/coin /usr/local/lib/coin/ThirdParty)
	
   FIND_LIBRARY(COIN_LAPACK_LIBRARY	  coinlapack HINTS  ${IPOPT_HOME}/coin/lib/ ${IPOPT_HOME}/coin/lib/ThirdParty/ ${IPOPT_HOME}/lib/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ /usr/local/lib/coin /usr/local/lib/coin/ThirdParty)
	
   FIND_LIBRARY(COIN_BLAS_LIBRARY  coinblas HINTS  ${IPOPT_HOME}/coin/lib/ ${IPOPT_HOME}/coin/lib/ThirdParty/ ${IPOPT_HOME}/lib/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ ${IPOPT_HOME}/lib/coin/ ${IPOPT_HOME}/lib/coin/ThirdParty/ /usr/local/lib/coin /usr/local/lib/coin/ThirdParty)
	
	  FIND_LIBRARY(MUMPS_LIBRARY HINTS ${IPOPT_HOME}/)
	
	  MESSAGE(STATUS "${PTHREAD_LIBRARIES}")
	
	  IF (IPOPT_LIBRARY)
	    SET(IPOPT_FOUND_LIBS TRUE)
	    SET(IPOPT_LIBRARIES ${IPOPT_LIBRARY})
	
	    IF(COIN_HSL_LIBRARY)
	      SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${COIN_HSL_LIBRARY})
	    ENDIF(COIN_HSL_LIBRARY)
	
	    IF(COIN_MUMPS_LIBRARY)
	      SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${COIN_MUMPS_LIBRARY})
	    ENDIF(COIN_MUMPS_LIBRARY)
	   
	    IF(COIN_METIS_LIBRARY)
	      SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${COIN_METIS_LIBRARY})
	    ENDIF(COIN_METIS_LIBRARY)
	
	    IF(COIN_LAPACK_LIBRARY AND COIN_BLAS_LIBRARY)
	      SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${COIN_LAPACK_LIBRARY} ${COIN_BLAS_LIBRARY})
	    ELSE(COIN_LAPACK_LIBRARY AND COIN_BLAS_LIBRARY)
	      SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${LAPACK_LIBRARIES})
	    ENDIF(COIN_LAPACK_LIBRARY AND COIN_BLAS_LIBRARY)
	         
	    SET(IPOPT_LIBRARIES ${IPOPT_LIBRARIES} ${EXTRA_LIBRARIES} ${CMAKE_DL_LIBS} ${CMAKE_THREAD_LIBS_INIT} -lgfortran)
	
	    MESSAGE(STATUS "Found Ipopt libs: ${IPOPT_LIBRARIES}")
	  ELSE (IPOPT_LIBRARY)
	    MESSAGE(STATUS "Could not find Ipopt libs")
	  ENDIF (IPOPT_LIBRARY)
	
	  IF(IPOPT_FOUND_INCLUDE AND IPOPT_FOUND_LIBS)
	    SET(IPOPT_FOUND TRUE)
	  ENDIF(IPOPT_FOUND_INCLUDE AND IPOPT_FOUND_LIBS)

   include(FindPackageHandleStandardArgs)
   # handle the QUIETLY and REQUIRED arguments and set LIBIPOPT_FOUND to TRUE
   # if all listed variables are TRUE
   find_package_handle_standard_args(IPOPT  DEFAULT_MSG
                                     IPOPT_LIBRARY IPOPT_INCLUDE_DIR)

   mark_as_advanced(IPOPT_INCLUDE_DIR IPOPT_LIBRARY )
   
ENDIF()