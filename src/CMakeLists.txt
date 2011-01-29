file(GLOB_RECURSE ESPRESSO_SOURCES *.cpp)
file(GLOB_RECURSE NOT_ESPRESSO_SOURCES
  "main/_espresso.cpp"
  "main/pypresso.cpp"
  "esutil/GenLogger.cpp"
  iterator/unittest/*.cpp
  interaction/unittest/*.cpp
  unittest/*.cpp
  analysis/unittest/*.cpp
  esutil/unittest/*.cpp
  storage/unittest/*.cpp
  bc/unittest/*.cpp
  integrator/unittest/*.cpp
)

list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.hpp.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/acconfig.hpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if (EXTERNAL_BOOST)
  set (BOOST_LIB ${Boost_LIBRARIES})
else(EXTERNAL_BOOST)
  set (BOOST espresso_boost)
endif(EXTERNAL_BOOST)

add_library(espresso_common ${ESPRESSO_SOURCES})
target_link_libraries(espresso_common ${BOOST} mpi4py ${PYTHON_LIBRARIES} ${MPI_LIBRARIES})

add_executable(pypresso main/pypresso.cpp)
target_link_libraries(pypresso espresso_common)

add_library(_espresso main/_espresso.cpp)
target_link_libraries(_espresso espresso_common)
