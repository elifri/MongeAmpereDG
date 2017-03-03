if (METHOD_OF_ELLIPSOIDS_INCLUDES)
  set(METHOD_OF_ELLIPSOIDS_QUIETLY TRUE)
endif (METHOD_OF_ELLIPSOIDS_INCLUDES)

find_path(METHOD_OF_ELLIPSOIDS_INCLUDES
  NAMES
  SupportingEllipsoids/MethodOfEllipsoids.hpp
  PATHS
  $ENV{IGPM_T2_LIBDIR}
  ${INCLUDE_INSTALL_DIR}
  $ENV{HOME}/workspace/Method\ of\ Ellipsoids/src
  ${METHOD_OF_ELLIPSOIDS_PATH}/src
)

message("searching for ${HOME}/workspace/Method\ of\ Ellipsoids/src/SupportingEllipsoids/MethodOfEllipsoids.hpp")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METHOD_OF_ELLIPSOIDS DEFAULT_MSG
                                  METHOD_OF_ELLIPSOIDS_INCLUDES)

mark_as_advanced(METHOD_OF_ELLIPSOIDS_INCLUDES)
