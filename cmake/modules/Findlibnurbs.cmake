if (LIBNURBS_INCLUDES)
  set(LIBNURBS_QUIETLY TRUE)
endif (LIBNURBS_INCLUDES)

find_path(LIBNURBS_INCLUDES
  NAMES
  opennurbs.h
  PATHS
  $ENV{LIBNURBS_DIR}
  ${INCLUDE_INSTALL_DIR}
  /home/disk/friebel/workspace/opennurbs_20130711
)

FIND_LIBRARY(LIBNURBS_LIBRARIES
  NAMES
  libopenNURBS.a
  PATHS
  $ENV{LIBNURBS_DIR}
  ${INCLUDE_INSTALL_DIR}
  /home/disk/friebel/workspace/opennurbs_20130711
)


message("searching for libnurbs .. in  $ENV{HOME}/workspace/opennurbs_20130711/opennurbs.h")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LIBNURBS DEFAULT_MSG
                                  LIBNURBS_INCLUDES LIBNURBS_LIBRARIES)

mark_as_advanced(LIBNURBS_INCLUDES LIBNURBS_LIBRARIES)
