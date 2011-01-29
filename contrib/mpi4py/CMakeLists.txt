set(MPI4PY_DIR mpi4py-1.2.1) 
set(MPI4PY_SOURCES MPI.cpp dynload.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/${MPI4PY_DIR}/src)

add_library(mpi4py ${MPI4PY_SOURCES})
target_link_libraries(mpi4py ${PYTHON_LIBRARIES} ${MPI_LIBRARIES})
