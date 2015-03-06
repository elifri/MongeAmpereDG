
if (IGPM_T2_LIB_INCLUDES)
  set(IGPM_T2_LIB_FIND_QUIETLY TRUE)
endif (IGPM_T2_LIB_INCLUDES)

find_path(IGPM_T2_LIB_INCLUDES
  NAMES
  include/igpm_t2_lib.hpp
  PATHS
  $ENV{IGPM_T2_LIBDIR}
  ${INCLUDE_INSTALL_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IGPM_T2_LIB DEFAULT_MSG
                                  IGPM_T2_LIB_INCLUDES)

mark_as_advanced(IGPM_T2_LIB_INCLUDES)
