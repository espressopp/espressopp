# - Find FFTW3
# Find the native FFTW3 includes and library
#
# FFTW3_INCLUDE_DIRS - where to find fftw3.h
# FFTW3_LIBRARIES - List of libraries when using FFTW3.
# FFTW3_FOUND - True if FFTW3 found.

if(FFTW3_INCLUDE_DIR AND FFTW3_LIBRARY)
  # Already in cache, be silent
  set (FFTW3_FIND_QUIETLY TRUE)
endif()

find_package(PkgConfig)
pkg_check_modules(PC_FFTW3 fftw3)
find_path(FFTW3_INCLUDE_DIR fftw3.h HINTS ${PC_FFTW3_INCLUDE_DIRS})

find_library(FFTW3_LIBRARY NAMES fftw3 HINTS ${PC_FFTW3_LIBRARY_DIRS} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE
# if all listed variables are TRUE

find_package_handle_standard_args(FFTW3 DEFAULT_MSG FFTW3_LIBRARY FFTW3_INCLUDE_DIR )

# Copy the results to the output variables and target.
if(FFTW3_FOUND)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARY} )
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR} )

  if(NOT TARGET FFTW3::fftw3)
    add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
    set_target_properties(FFTW3::fftw3 PROPERTIES
      IMPORTED_LINK_INTERFACE_LANGUAGES "C"
      IMPORTED_LOCATION "${FFTW3_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}")
  endif()
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY )
