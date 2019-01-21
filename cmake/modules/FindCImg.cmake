
if (CIMG_INCLUDES)
  set(CIMG_QUIETLY TRUE)
endif (CIMG_INCLUDES)

IF (CIMG_FOUND)
ELSE (CIMG_FOUND)
  find_path(CIMG_INCLUDES
    NAMES
    CImg.h
    PATHS
    ${INCLUDE_INSTALL_DIR}
    $ENV{HOME}/workspace/CImg-1.6.0/
    ${CIMG_PATH}
    ../CImg-1.6.0
  )
  
  message("searching for $ENV{HOME}/workspace/Cimg-1.6.0/CImg.hpp and ${CIMG_PATH}/CImg.hpp"  )
  
  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(SCIMG DEFAULT_MSG
                                    CIMG_INCLUDES)
  
  mark_as_advanced(CIMG_INCLUDES)
  
  IF( CIMG_FOUND AND NOT CIMG_QUIETLY )
      MESSAGE(STATUS "package CImg found")
  ENDIF()
ENDIF(CIMG_FOUND)
