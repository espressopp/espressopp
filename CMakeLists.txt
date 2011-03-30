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
option(BUILD_STANDALONE "Builds pypresso an standalone binary" OFF)
option(BUILD_PYTHON_MODULE "Builds epresso as python module" ON)

if(NOT BUILD_STANDALONE AND NOT BUILD_PYTHON_MODULE)
  message(FATAL_ERROR "You need to enable pypresso or the python module")
endif(NOT BUILD_STANDALONE AND NOT BUILD_PYTHON_MODULE)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

########################################################################
#Process MPI settings
########################################################################

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

########################################################################
#Process Python settings
########################################################################

find_package(PythonInterp REQUIRED)

execute_process(COMMAND ${PYTHON_EXECUTABLE}
                -c "import distutils.sysconfig as cg; print cg.get_python_inc()"
                OUTPUT_VARIABLE PYTHON_INCLUDE_PATH
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set(PYTHON_INCLUDE_PATH ${PYTHON_INCLUDE_PATH} CACHE PATH "Python Include Directory")
mark_as_advanced(PYTHON_INCLUDE_PATH)

message(STATUS "PYTHON_INCLUDE_PATH = ${PYTHON_INCLUDE_PATH}")
include_directories(${PYTHON_INCLUDE_PATH})

execute_process(COMMAND ${PYTHON_EXECUTABLE} 
                -c "import distutils.sysconfig as cg; print cg.get_python_lib(1,0,prefix='${CMAKE_INSTALL_EXEC_PREFIX}')"
                OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${PYTHON_EXECUTABLE}
                        -c "import sys; print sys.version[:3]"
                        OUTPUT_VARIABLE PYTHON_VERSION
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "PYTHON_VERSION = ${PYTHON_VERSION}")

execute_process(COMMAND ${PYTHON_EXECUTABLE}
                -c "import distutils.sysconfig as cg; print cg.get_config_var('LIBDIR')"
                OUTPUT_VARIABLE PYTHON_LIBDIR
                OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "PYTHON_LIBDIR = ${PYTHON_LIBDIR}")

find_library(PYTHON_LIBRARIES
             NAMES python${PYTHON_VERSION}
             PATHS ${PYTHON_LIBDIR}  )

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")

mark_as_advanced(PYTHON_INCLUDE_PATH PYTHON_LIBRARIES)

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
