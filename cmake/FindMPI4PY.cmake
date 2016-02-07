# - Find MPI4PY
# Find the native MPI4PY includes and library
#
# MPI4PY_INCLUDES - where to find mpi4py.h
# MPI4PY_LIBRARIES - List of libraries when using MPI4PY.
# MPI4PY_FOUND - True if MPI4PY found.

if (MPI4PY_INCLUDES)
  # Already in cache, be silent
  set (MPI4PY_FIND_QUIETLY TRUE)
endif (MPI4PY_INCLUDES)

if(PYTHON_EXECUTABLE)
  execute_process(COMMAND ${PYTHON_EXECUTABLE} 
                -c "import distutils.sysconfig as cg; print cg.get_python_lib(1,0)"
		OUTPUT_VARIABLE PYTHON_SITEDIR OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND "${PYTHON_EXECUTABLE}" 
                -c "import sys, mpi4py; sys.stdout.write(mpi4py.__version__)"
                    OUTPUT_VARIABLE MPI4PY_VERSION
                    RESULT_VARIABLE _MPI4PY_VERSION_RESULT
                    ERROR_QUIET)
  if(NOT _MPI4PY_VERSION_RESULT)
    message("-- mpi4py version: " ${MPI4PY_VERSION})
  else()
    set(MPI4PY_VERSION 0.0)
  endif()

  execute_process(COMMAND
      "${PYTHON_EXECUTABLE}" "-c" "import mpi4py; print mpi4py.get_include()"
      OUTPUT_VARIABLE MPI4PY_INCLUDE_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  
endif(PYTHON_EXECUTABLE)

find_path (MPI4PY_INCLUDES mpi4py/mpi4py.h HINTS ${MPI4PY_INCLUDE_DIR} ${PYTHON_SITEDIR}/mpi4py/include )
if(NOT MPI4PY_INCLUDES)
  message("     mpi4py.h not found. Please make sure you have installed the developer version of mpi4py")
endif()

find_file (MPI4PY_LIBRARIES MPI.so HINTS ${MPI4PY_INCLUDE_DIR}/.. ${PYTHON_SITEDIR}/mpi4py)

# handle the QUIETLY and REQUIRED arguments and set MPI4PY_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPI4PY REQUIRED_VARS MPI4PY_LIBRARIES MPI4PY_INCLUDES VERSION_VAR MPI4PY_VERSION)

mark_as_advanced (MPI4PY_LIBRARIES MPI4PY_INCLUDES)
