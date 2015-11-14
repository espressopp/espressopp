if (NOT EXTERNAL_BOOST)
  add_subdirectory(boost)
endif(NOT EXTERNAL_BOOST)
if(NOT EXTERNAL_MPI4PY)
add_subdirectory(mpi4py)
endif(NOT EXTERNAL_MPI4PY)

