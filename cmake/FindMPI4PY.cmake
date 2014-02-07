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
endif(PYTHON_EXECUTABLE)

find_path (MPI4PY_INCLUDES mpi4py/mpi4py.h HINTS ${PYTHON_SITEDIR}/mpi4py/include )

find_file (MPI4PY_LIBRARIES MPI.so HINTS ${PYTHON_SITEDIR}/mpi4py)

# handle the QUIETLY and REQUIRED arguments and set MPI4PY_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (MPI4PY DEFAULT_MSG MPI4PY_LIBRARIES MPI4PY_INCLUDES)

mark_as_advanced (MPI4PY_LIBRARIES MPI4PY_INCLUDES)
