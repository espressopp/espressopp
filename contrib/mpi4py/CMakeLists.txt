enable_language(C)

set(MPI4PY_DIR mpi4py-3.0.3) 
set(MPI4PY_SOURCES ${MPI4PY_DIR}/src/MPI.c)

add_library(MPI ${MPI4PY_SOURCES})
target_link_libraries(MPI PUBLIC Python3::Python MPI::MPI_CXX)
target_include_directories(MPI PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${MPI4PY_DIR}/src>)
set_target_properties(MPI PROPERTIES PREFIX "" SUFFIX ".so")

install(TARGETS MPI LIBRARY DESTINATION ${PYTHON_INSTDIR}/mpi4py ARCHIVE DESTINATION ${PYTHON_INSTDIR}/mpi4py)
install(FILES __init__.py DESTINATION ${PYTHON_INSTDIR}/mpi4py)
if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_command(TARGET MPI COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py)
endif()

add_library(MPI4PY_include INTERFACE)
target_include_directories(MPI4PY_include SYSTEM INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${MPI4PY_DIR}/src/mpi4py/include>)

add_library(MPI4PY::mpi4py ALIAS MPI)
add_library(MPI4PY::include ALIAS MPI4PY_include)
