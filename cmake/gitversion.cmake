######################################
# Create gitversion.hpp file         #
######################################
set(CMAKE_MODULE_PATH ${TOP_SOURCE_DIR}/cmake)
find_package(Git)
if (GIT_FOUND AND IS_DIRECTORY ${TOP_SOURCE_DIR}/.git)
  execute_process( COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY ${TOP_SOURCE_DIR}
    OUTPUT_VARIABLE ESPP_GIT_ID OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process( COMMAND ${GIT_EXECUTABLE} diff-index --name-only HEAD
    WORKING_DIRECTORY ${TOP_SOURCE_DIR}
    OUTPUT_VARIABLE _HAS_CHANGES OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if (NOT "${_HAS_CHANGES}" STREQUAL "")
    set(ESPP_GIT_ID "${ESPP_GIT_ID} (dirty)")
  endif()
  message("Current revision is ${ESPP_GIT_ID}")
  set (ESPP_GIT_ID ${ESPP_GIT_ID})
else (GIT_FOUND AND IS_DIRECTORY ${TOP_SOURCE_DIR}/.git)
  set (ESPP_GIT_ID)
endif (GIT_FOUND AND IS_DIRECTORY ${TOP_SOURCE_DIR}/.git)
set (GIT_HEADER "gitversion.hpp")
set (NEW_GIT_HEADER "new_gitversion.hpp")
file(WRITE ${NEW_GIT_HEADER} "static const std::string gitversion = \"${ESPP_GIT_ID}\";\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${NEW_GIT_HEADER} ${GIT_HEADER})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${NEW_GIT_HEADER})

