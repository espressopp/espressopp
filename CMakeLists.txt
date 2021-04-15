cmake_minimum_required(VERSION 3.12)

project(espressopp LANGUAGES CXX)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

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
option(WITH_XTC "Build with DumpXTC class (requires libgromacs)" OFF)
option(BUILD_SHARED_LIBS "Build shared libs" ON)
if(NOT BUILD_SHARED_LIBS)
  message(WARNING "Building static libraries might lead to problems with python modules - you are on your own!")
endif()

option(USE_GCOV "Enable gcov support" OFF)
if(USE_GCOV)
  message(STATUS "Enabling gcov support")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} --coverage -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage -O0")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  SET(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

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

find_package(MPI REQUIRED COMPONENTS CXX)

########################################################################
#Process FFTW3 settings
########################################################################

find_package(FFTW3 REQUIRED)
include_directories(${FFTW3_INCLUDES})

########################################################################
#Process GROMACS settings
########################################################################

if(WITH_XTC)
  find_package(GROMACS REQUIRED)
  include_directories(${GROMACS_INCLUDE_DIRS})
  add_definitions( -DHAS_GROMACS ${GROMACS_DEFINITIONS})
endif()


########################################################################
#Process Python settings
########################################################################

find_package(Python3 COMPONENTS Interpreter Development)

execute_process(COMMAND ${Python3_EXECUTABLE}
                -c "import distutils.sysconfig as cg; print(cg.get_python_inc())"
                OUTPUT_VARIABLE PYTHON_INCLUDE_PATH
                OUTPUT_STRIP_TRAILING_WHITESPACE)

set(PYTHON_INCLUDE_PATH ${Python3_INCLUDE_DIRS} CACHE PATH "Python Include Directory")
mark_as_advanced(PYTHON_INCLUDE_PATH)

message(STATUS "PYTHON_INCLUDE_PATH = ${PYTHON_INCLUDE_PATH}")
include_directories(${PYTHON_INCLUDE_PATH})

if(NOT PYTHON_INSTDIR)
    execute_process(COMMAND ${Python3_EXECUTABLE}
                -c "import distutils.sysconfig as cg; print(cg.get_python_lib(1,0,prefix='${CMAKE_INSTALL_EXEC_PREFIX}'))"
                OUTPUT_VARIABLE PYTHON_INSTDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

execute_process(COMMAND ${Python3_EXECUTABLE}
                        -c "import sys; print(sys.version[:3])"
                        OUTPUT_VARIABLE PYTHON_VERSION
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(COMMAND ${Python3_EXECUTABLE}
                        -c "import sys; print(sys.version[:3].replace('.', ''))"
                        OUTPUT_VARIABLE PYTHON_VERSION_NO_DOT
                        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "PYTHON_VERSION = ${PYTHON_VERSION}")
message(STATUS "PYTHON_LIBDIR = ${Python3_LIBRARY_DIRS}")

set(PYTHON_LIBRARIES ${Python3_LIBRARIES})

if(NOT Python3_LIBRARIES)
  message(FATAL_ERROR "Python libraries not found!")
endif(NOT Python3_LIBRARIES)
message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")

mark_as_advanced(PYTHON_INCLUDE_PATH PYTHON_LIBRARIES)

########################################################################
# Basic system tests (standard libraries, headers, functions, types)   #
########################################################################
list(APPEND EXTRA_LIBRARIES ${CMAKE_DL_LIBS})

# Minimum required Boost version is set to 1.69.0
set(NEEDED_Boost_VERSION 1.69.0)
find_package(Boost ${NEEDED_Boost_VERSION} REQUIRED COMPONENTS mpi serialization system filesystem python${PYTHON_VERSION_NO_DOT})
set(Boost_LIBS ${Boost_LIBRARIES})

message(STATUS "BOOST_LIBS = ${Boost_LIBRARIES}")

# Append Python library to the list of Boost libraries.
list(APPEND Boost_LIBS ${Boost_LIBRARIES})
set(Boost_LIBRARIES ${Boost_LIBS})

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost (required at least version ${NEEDED_Boost_VERSION}) not found, make sure you have installed boost and it's dev packages (mpi, serialization, system, filesystem, python3). To hint cmake to your Boost installation, use the option -DBOOST_ROOT=/path/to/boost.")
endif(NOT Boost_FOUND)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${CMAKE_SOURCE_DIR}/src)

option(EXTERNAL_MPI4PY "Use external mpi4py" ON)

if(EXTERNAL_MPI4PY)
  find_package(MPI4PY 3.0.0 REQUIRED)
  include_directories(${MPI4PY_INCLUDES})
else(EXTERNAL_MPI4PY)
  include_directories(${CMAKE_SOURCE_DIR}/contrib/mpi4py/mpi4py-3.0.3/src/mpi4py/include)
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
add_subdirectory(contrib)
set (TEST_ENV "PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")
add_subdirectory(testsuite)

add_custom_target(symlink ALL COMMENT "Creating symlink")
add_custom_command(TARGET symlink COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_SOURCE_DIR}/src espressopp)

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

######################################
# User documentation (Sphinx)     #
######################################

find_program(SPHINX_EXECUTABLE
  NAMES sphinx-build
  DOC "Sphinx documentation generation tool (http://www.sphinx-doc.org/)"
)

if ( SPHINX_EXECUTABLE )

   set ( SPHINX_BUILD_ROOT "${CMAKE_SOURCE_DIR}/doc/ug" )

   add_custom_target(ug
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV}" "${SPHINX_EXECUTABLE}" -b html . _build/html
      COMMENT "HTML ug documentation available at ${SPHINX_BUILD_ROOT}/_build/html"
      DEPENDS ${SPHINX_BUILD_ROOT}/conf.py _espressopp
      WORKING_DIRECTORY ${SPHINX_BUILD_ROOT}
   )

   message ( STATUS "Generate user documentation via ${SPHINX_EXECUTABLE}: make ug")

   add_custom_target(ug-pdf
      COMMAND ${CMAKE_COMMAND} -E env "${TEST_ENV}" "${SPHINX_EXECUTABLE}" -b latex . _build/latex
      COMMAND make -C "${SPHINX_BUILD_ROOT}/_build/latex" all-pdf
      COMMENT "TeX ug documentation available at ${SPHINX_BUILD_ROOT}/_build/latex"
      DEPENDS ${SPHINX_BUILD_ROOT}/conf.py _espressopp
      WORKING_DIRECTORY ${SPHINX_BUILD_ROOT}
   )

   message ( STATUS "Generate LaTeX manual via ${SPHINX_EXECUTABLE}: make ug-pdf")
else ( SPHINX_EXECUTABLE )
    message ( WARNING "Not building user documentation because Sphinx not found.")
endif ( SPHINX_EXECUTABLE )
