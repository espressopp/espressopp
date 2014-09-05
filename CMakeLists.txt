cmake_minimum_required(VERSION 2.6)

project(espressopp)

set(PROJECT_VERSION "0.9.0-dev")

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
option(EXTERNAL_BOOST "Use external boost" ON)
option(BUILD_SHARED_LIBS "Build shared libs" ON)
option(BUILD_STANDALONE "Builds pypresso an standalone binary" OFF)
option(BUILD_PYTHON_MODULE "Builds epresso as python module" ON)
if (NOT DEFINED LIB)
  set(LIB "lib")
endif(NOT DEFINED LIB)

if(NOT BUILD_STANDALONE AND NOT BUILD_PYTHON_MODULE)
  message(FATAL_ERROR "You need to enable pypresso or the python module")
endif(NOT BUILD_STANDALONE AND NOT BUILD_PYTHON_MODULE)

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

########################################################################
#VampirTrace settings
########################################################################

set(VAMPIRTRACE False)

if(VAMPIRTRACE)
   set(VAMPIRTRACE_ROOT /people/thnfs/homes/stuehn/software)
   set(VAMPIRTRACE_LIBRARIES ${VAMPIRTRACE_ROOT}/lib/libvt-mpi.so
                             ${VAMPIRTRACE_ROOT}/lib/libotf.so)
   # Attention: VampirTrace should be included before MPI
   include_directories(${VAMPIRTRACE_ROOT}/include)
   add_definitions( -DVTRACE )
endif(VAMPIRTRACE)

########################################################################
#Process MPI settings
########################################################################

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

########################################################################
#Process FFTW3 settings
########################################################################

find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDES})

########################################################################
#Process Python settings
########################################################################

find_package(PythonInterp 2 REQUIRED)

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

# this conflicts with mpi4py-1.3 (is used to check for MPI version there)
# add_definitions( -DHAVE_CONFIG_H )

if (EXTERNAL_BOOST)
  find_package(Boost 1.45.0 COMPONENTS mpi python serialization filesystem)
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost not found, make sure you have installed boost and it's dev packages. To hint cmake to your Boost installation, use the option -DBOOST_ROOT=/path/to/boost. As last resort you can use our internal replacement with -DEXTERNAL_BOOST=OFF")
  endif(NOT Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else(EXTERNAL_BOOST)
  include_directories(${CMAKE_SOURCE_DIR}/contrib/boost)
endif(EXTERNAL_BOOST)
include_directories(${CMAKE_SOURCE_DIR}/src)

option(EXTERNAL_MPI4PY "Use external mpi4py" ON)

if(EXTERNAL_MPI4PY)
  find_package(MPI4PY 1.3 REQUIRED)
  include_directories(${MPI4PY_INCLUDES})
else(EXTERNAL_MPI4PY)
  include_directories(${CMAKE_SOURCE_DIR}/contrib/mpi4py/mpi4py-1.3/src/include)
  set(MPI4PY_LIBRARIES MPI)
endif(EXTERNAL_MPI4PY)

# Required for acconfig.h
include_directories(${CMAKE_BINARY_DIR}/src)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set ( RT_LIBRARIES )
else()
  set ( RT_LIBRARIES rt )
endif()

######################################
# Include the following subdirectory # 
######################################
add_subdirectory(src)
add_subdirectory(contrib)

add_custom_target(symlink ALL COMMENT "Creating symlink")
add_custom_command(TARGET symlink COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_BINARY_DIR}/src espresso)

######################################
# System documentation (Doxygen)     # 
######################################

find_package( Doxygen QUIET )

if ( DOXYGEN_FOUND )

    set ( DOXYGEN_BUILD_ROOT "${CMAKE_SOURCE_DIR}/doc/dg" )
    file ( MAKE_DIRECTORY ${DOXYGEN_BUILD_ROOT} )

   # The initial rm command gets rid of everything previously built by this
   # custom command.

   add_custom_command(
      OUTPUT ${DOXYGEN_BUILD_ROOT}/html/index.html
      COMMAND rm -rf ${DOXYGEN_BUILD_ROOT}/html
      COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_BUILD_ROOT}/doxygen.cfg
      COMMAND echo "HTML documentation available at ${DOXYGEN_BUILD_ROOT}/html"
      DEPENDS ${DOXYGEN_BUILD_ROOT}/doxygen.cfg
      WORKING_DIRECTORY ${DOXYGEN_BUILD_ROOT}
   )

   add_custom_target(
      doc
      DEPENDS
      ${DOXYGEN_BUILD_ROOT}/html/index.html
   )

   message ( STATUS "Generate system documentation via ${DOXYGEN_EXECUTABLE}: make doc")

else ( DOXYGEN_FOUND )
    message ( WARNING "Not building system documentation because Doxygen not found.")
endif ( DOXYGEN_FOUND )

