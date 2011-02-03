cmake_minimum_required(VERSION 2.6)

project(ESPResSo++)

set(PROJECT_VERSION "0.9")

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

########################################################################
#Process MPI settings
########################################################################
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_PATH})

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
include(CheckIncludeFiles)
check_include_files(sys/time.h   HAVE_SYS_TIME_H)
check_include_files(sys/resource.h HAVE_SYS_RESOURCE_H)
check_include_files(unistd.h     HAVE_UNISTD_H)

list(APPEND EXTRA_LIBRARIES ${CMAKE_DL_LIBS})

add_definitions( -DHAVE_CONFIG_H )

if (EXTERNAL_BOOST)
  find_package(Boost REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
  target_link_libraries(foo ${Boost_LIBRARIES})
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

