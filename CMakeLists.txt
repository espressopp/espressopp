cmake_minimum_required(VERSION 2.6)

project(espressopp)

set(PROJECT_VERSION "0.9.0-dev")
SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "9")
SET(CPACK_PACKAGE_VERSION_PATCH "0")

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


enable_language(C)
enable_language(CXX)


######################################################################
# cpack seetings
#####################################################################
set(CPACK_CMAKE_GENERATOR "Unix Makefiles")
set(CPACK_GENERATOR "")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_VENDOR "espresso-pp.de")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PROJECT_NAME}-${PROJECT_VERSION}")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "ESPResSo++ is a highly versatile software package for the scientific simulation ")
set(CPACK_SOURCE_IGNORE_FILES ".*~;/\\\\.hg/;\\\\.hgignore;/CMakeFiles/;CMakeCache\\\\.txt$;cmake_install\\\\.cmake$;CPack.*\\\\.cmake$;Makefile$;ESPRC$")
if (EXTERNAL_BOOST)
  set(CPACK_SOURCE_IGNORE_FILES "${CPACK_SOURCE_IGNORE_FILES};contrib/boost/")
endif (EXTERNAL_BOOST)
include(CPack)


######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

########################################################################
# User input options                                                   #
########################################################################
option(EXTERNAL_BOOST "Use external boost" OFF)
option(BUILD_SHARED_LIBS "Build shared libs" ON)
option(BUID_STANDALONE "Builds pypresso an standalone binary" OFF)
option(BUID_PYTHON_MODULE "Builds epresso as python module" ON)

if(NOT BUID_STANDALONE AND NOT BUID_PYTHON_MODULE)
  message(FATAL_ERROR "You need to enable pypresso or the python module")
endif(NOT BUID_STANDALONE AND NOT BUID_PYTHON_MODULE)


########################################################################
#Process MPI settings
########################################################################
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIR})

find_package(PythonInterp REQUIRED)
execute_process(COMMAND ${PYTHON_EXECUTABLE} -c "from distutils import sysconfig;print sysconfig.get_python_lib(1,0,prefix='${CMAKE_INSTALL_EXEC_PREFIX}')"
    OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFiles)
check_include_files("sys/time.h;sys/resource.h;unistd.h" CMAKE_HEADERS)

list(APPEND EXTRA_LIBRARIES ${CMAKE_DL_LIBS})

add_definitions( -DHAVE_CONFIG_H )

if (EXTERNAL_BOOST)
  find_package(Boost 1.42.0 REQUIRED COMPONENTS mpi python serialization )
  include_directories(${Boost_INCLUDE_DIRS})
else(EXTERNAL_BOOST)
  include_directories(${CMAKE_SOURCE_DIR}/contrib/boost)
endif(EXTERNAL_BOOST)
include_directories(${CMAKE_SOURCE_DIR}/src)

# Required for acconfig.h
include_directories(${CMAKE_BINARY_DIR}/src)

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(src)
add_subdirectory(contrib)

add_custom_target(symlink ALL COMMENT "Creating symlink")
add_custom_command(TARGET symlink COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_BINARY_DIR}/src espresso)
