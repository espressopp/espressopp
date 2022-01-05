cmake_minimum_required(VERSION 3.12)

include(CMakePrintHelpers)

project(espressopp LANGUAGES C CXX)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(ESPP_CMAKE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${ESPP_CMAKE_DIR})

enable_testing()

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCXXCompilerFlag)

# Check for C++14 standard
set(CMAKE_CXX_STANDARD 14) # C++14...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

########################################################################
# User input options                                                   #
########################################################################
option(ESPP_LOCAL_ARCHITECTURE "Use instruction set of the local architecture." OFF)
option(ESPP_VEC_REPORT "Enable reporting of loop vectorization." OFF)
option(ESPP_WERROR "Treat warnings as errors." OFF)
option(ESPP_WALL "Build with more warnings." ON)
option(BUILD_SHARED_LIBS "Build shared libs" ON)
if(NOT BUILD_SHARED_LIBS)
    message(WARNING "Building static libraries might lead to problems with python modules - you are on your own!")
endif()

option(USE_GCOV "Enable gcov support" OFF)
if(USE_GCOV)
    message(STATUS "Enabling gcov support")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif()

########################################################################
#Process MPI settings
########################################################################

find_package(MPI REQUIRED COMPONENTS CXX)

########################################################################
#Process FFTW3 settings
########################################################################

find_package(FFTW3 REQUIRED)

########################################################################
#Find HDF5
########################################################################

set(HDF5_PREFER_PARALLEL ON)
#set(HDF5_FIND_DEBUG ON)
find_package(HDF5 REQUIRED COMPONENTS C HL)
if (NOT HDF5_IS_PARALLEL)
    message(FATAL_ERROR "HDF5 with MPI support is required! Serial version detected.")
endif()

########################################################################
#Process GROMACS settings
########################################################################
find_package(GROMACS 2021 CONFIG QUIET NAMES gromacs gromacs_d)
option(WITH_XTC "Build with DumpXTC class (requires libgromacs)" ${GROMACS_FOUND})
if(WITH_XTC)
	find_package(GROMACS 2021 CONFIG REQUIRED NAMES gromacs gromacs_d)
endif()


########################################################################
#Process Python settings
########################################################################

find_package(Python3 COMPONENTS Interpreter Development)
#cmake_print_variables( Python3_INCLUDE_DIRS Python3_LINK_OPTIONS Python3_LIBRARY_DIRS Python3_RUNTIME_LIBRARY_DIRS Python3_VERSION Python3_VERSION_MAJOR Python3_VERSION_MINOR Python3_NumPy_VERSION)
set(PYTHON_VERSION_NO_DOT "${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}")

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
list(APPEND EXTRA_LIBRARIES ${CMAKE_DL_LIBS})

#set( Boost_DEBUG ON )
find_package(Boost 1.69.0 REQUIRED COMPONENTS mpi serialization system filesystem python${PYTHON_VERSION_NO_DOT} numpy${PYTHON_VERSION_NO_DOT} unit_test_framework)
#cmake_print_variables( Boost_INCLUDE_DIRS Boost_LIBRARY_DIRS Boost_LIBRARIES Boost_LIB_VERSION  )

######################################
# Automatic CCache detection         #
######################################
option(USE_CCACHE "automatically search for and use ccache" OFF)
if (USE_CCACHE)
    find_program(CCACHE_PROGRAM ccache)
    if(CCACHE_PROGRAM)
        message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
        set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    endif()
endif()

######################################
# Include the following subdirectory #
######################################
add_subdirectory(src)
add_subdirectory(testsuite)

add_custom_target(symlink ALL COMMENT "Creating symlink")
add_custom_command(TARGET symlink COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_SOURCE_DIR}/src espressopp)

add_subdirectory(doc)
