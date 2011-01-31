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

file(GLOB_RECURSE PYTHON_SCRIPTS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)

if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(scripts ALL COMMENT "Copying scripts to builddir: *.py")
  foreach(FILENAME ${PYTHON_SCRIPTS})
    set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
    set(DST "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
    add_custom_command(TARGET scripts COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST})
  endforeach(FILENAME)
endif(NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})

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
add_dependencies(pypresso scripts)

add_library(_espresso main/_espresso.cpp)
target_link_libraries(_espresso espresso_common)
add_dependencies(_espresso scripts)
