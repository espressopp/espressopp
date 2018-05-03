cmake_minimum_required(VERSION 3.1)

project(espressopp)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

enable_language(C)
enable_language(CXX)

enable_testing()

######################################################################
# compiler tests
# these need ot be done early (before further tests).
#####################################################################

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)

# Check for C++11 standard
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

########################################################################
# User input options                                                   #
########################################################################
option(EXTERNAL_BOOST "Use external boost" ON)
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

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

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
  add_definitions( -DHAS_GROMACS )
endif()


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

if(NOT PYTHON_LIBRARIES)
  message(FATAL_ERROR "Python libraries not found!")
endif(NOT PYTHON_LIBRARIES)
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

# Minimum required Boost version is set to 1.58 (only because of travis)
# but the version in contrib is 1.61;

if (EXTERNAL_BOOST)
  set(NEEDED_Boost_VERSION 1.58.0)
  find_package(Boost ${NEEDED_Boost_VERSION} REQUIRED COMPONENTS mpi python serialization system filesystem)
  if(NOT Boost_FOUND)
    message(FATAL_ERROR "Boost (required at least version ${NEEDED_Boost_VERSION}) not found, make sure you have installed boost and it's dev packages. To hint cmake to your Boost installation, use the option -DBOOST_ROOT=/path/to/boost. As last resort you can use our internal replacement with -DEXTERNAL_BOOST=OFF")
  endif(NOT Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
else(EXTERNAL_BOOST)
  include_directories(${CMAKE_SOURCE_DIR}/contrib/boost)
  set(Boost_LIBRARIES espressopp_boost)
endif(EXTERNAL_BOOST)
include_directories(${CMAKE_SOURCE_DIR}/src)

option(EXTERNAL_MPI4PY "Use external mpi4py" ON)

if(EXTERNAL_MPI4PY)
  find_package(MPI4PY 1.3.1 REQUIRED)
  include_directories(${MPI4PY_INCLUDES})
else(EXTERNAL_MPI4PY)
  include_directories(${CMAKE_SOURCE_DIR}/contrib/mpi4py/mpi4py-2.0.0/src/include)
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
if (EXTERNAL_MPI4PY)
  set (TEST_ENV "PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")
else (EXTERNAL_MPI4PY)
  set (TEST_ENV "PYTHONPATH=${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/contrib:$ENV{PYTHONPATH}")
endif (EXTERNAL_MPI4PY)
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
