# - Find DL
# Find the DL includes and library
#
# DL_INCLUDES - where to find dlfcn.h
# DL_LIBRARIES - List of libraries when using DL.
# DL_FOUND - True if DL found.

if (DL_INCLUDES)
  # Already in cache, be silent
  set (DL_FIND_QUIETLY TRUE)
endif (DL_INCLUDES)

find_path (DL_INCLUDES dlfcn.h)

find_library(DL_LIBRARIES NAMES dl)

# handle the QUIETLY and REQUIRED arguments and set DL_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args (DL DEFAULT_MSG DL_LIBRARIES DL_INCLUDES)

mark_as_advanced(DL_LIBRARIES DL_INCLUDES)
