# Random123
# - Check for the presence of RANDOM123
#
# The following variables are set when RANDOM123 is found:
#  RANDOM123_FOUND      = Set to true, if all components of RANDOM123 have been found.
#  RANDOM123_INCLUDES    = Include path for the header files of RANDOM123

find_path(RANDOM123_INCLUDE_DIR Random123/threefry.h)
if(RANDOM123_INCLUDE_DIR)
    set(RANDOM123_FOUND TRUE)
    message(STATUS "Found Random123: ${RANDOM123_INCLUDES}")
endif()

if (NOT RANDOM123_FOUND)

  if (NOT RANDOM123_ROOT_DIR)
    set (RANDOM123_ROOT_DIR ${CMAKE_INSTALL_PREFIX})
  endif (NOT RANDOM123_ROOT_DIR)

  find_path (RANDOM123_INCLUDES
    NAMES Random123/threefry.h Random123/u01.h
    HINTS ${RANDOM123_ROOT_DIR} ${CMAKE_INSTALL_PREFIX}
    PATH_SUFFIXES include
    )

  find_package_handle_standard_args (RANDOM123 DEFAULT_MSG RANDOM123_INCLUDES)

  if (RANDOM123_FOUND)
    if (NOT RANDOM123_FIND_QUIETLY)
      message (STATUS "Found components for RANDOM123")
      message (STATUS "RANDOM123_ROOT_DIR  = ${RANDOM123_ROOT_DIR}")
      message (STATUS "RANDOM123_INCLUDES  = ${RANDOM123_INCLUDES}")
    endif (NOT RANDOM123_FIND_QUIETLY)
    
    # Mark advanced variables
    mark_as_advanced (RANDOM123_ROOT_DIR RANDOM123_INCLUDES)
  else (RANDOM123_FOUND)
    if (RANDOM123_FIND_REQUIRED)
      message (FATAL_ERROR "Required package Random123 not found")
    endif (RANDOM123_FIND_REQUIRED)
  endif ()
endif()
