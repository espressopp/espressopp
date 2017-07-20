# - Find GROMACS
# Find the native GROMACS includes and library
#
# GROMACS_INCLUDES - where to find .h
# GROMACS_LIBRARIES - List of libraries when using GROMACS.
# GROMACS_FOUND - True if GROMACS found.

if (GROMACS_INCLUDES)
  # Already in cache, be silent
  set (GROMACS_FIND_QUIETLY TRUE)
endif (GROMACS_INCLUDES)

find_path (GROMACS_INCLUDES gromacs/fileio/xtcio.h)

find_library (GROMACS_LIBRARIES NAMES gromacs)

# handle the QUIETLY and REQUIRED arguments and set GROMACS_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (GROMACS DEFAULT_MSG GROMACS_LIBRARIES GROMACS_INCLUDES)

mark_as_advanced (GROMACS_LIBRARIES GROMACS_INCLUDES)
