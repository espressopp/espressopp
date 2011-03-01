file(GLOB_RECURSE ESPRESSO_SOURCES *.cpp)
file(GLOB_RECURSE NOT_ESPRESSO_SOURCES
  main/_espresso.cpp
  main/pypresso.cpp
  esutil/GenLogger.cpp
  esutil/PyLogger.cpp
  iterator/unittest/*.cpp
  interaction/unittest/*.cpp
  unittest/*.cpp
  analysis/unittest/*.cpp
  esutil/unittest/*.cpp
  storage/unittest/*.cpp
  bc/unittest/*.cpp
  integrator/unittest/*.cpp
  thermostat.cpp
)
list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})

file(GLOB_RECURSE PYTHON_SCRIPTS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)

add_custom_target(scripts ALL COMMENT "Copying scripts to builddir: *.py")
foreach(FILENAME ${PYTHON_SCRIPTS})
  set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
  set(DST "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
  string(REGEX MATCH .*/ DIR "${FILENAME}")
  if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_command(TARGET scripts COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST})
  endif(NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  install(FILES ${DST} DESTINATION ${PYTHON_INSTDIR}/espresso/${DIR})
endforeach(FILENAME)

#remove this two lines if autotools are gone
add_custom_command(TARGET scripts COMMAND ${CMAKE_COMMAND} -E touch 
  ${CMAKE_CURRENT_BINARY_DIR}/main/_setupPath.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/main/_setupPath.py DESTINATION ${PYTHON_INSTDIR}/espresso/main)

list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.hpp.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/acconfig.hpp)

#build ESPRC
#CMAKE_BINARY_DIR for espresso, CMAKE_CURRENT_BINARY_DIR for for _espressos.so and contrib for MPI.so
set (PYTHON_DIR "${CMAKE_BINARY_DIR}:${CMAKE_CURRENT_BINARY_DIR}/pymodule:${CMAKE_BINARY_DIR}/contrib/mpi4py")
#CMAKE_CURRENT_BINARY_DIR for pypresso
set (BIN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
#CMAKE_CURRENT_BINARY_DIR for libespresso_common.so and contib for libespresso_boost
set (LIB_DIR "${CMAKE_CURRENT_BINARY_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ESPRC.cmakein
  ${CMAKE_BINARY_DIR}/ESPRC @ONLY)

#after install ESPRC
set (PYTHON_DIR "${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR}")
set (BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set (LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ESPRC.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/ESPRC @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ESPRC DESTINATION bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if (EXTERNAL_BOOST)
  set (BOOST ${Boost_MPI_LIBRARY} ${Boost_PYTHON_LIBRARY} ${Boost_SERIALIZATION_LIBRARY} )
else(EXTERNAL_BOOST)
  set (BOOST espresso_boost)
endif(EXTERNAL_BOOST)

add_library(espresso_common ${ESPRESSO_SOURCES})
target_link_libraries(espresso_common ${BOOST} ${PYTHON_LIBRARIES} ${MPI_LIBRARIES})
install(TARGETS espresso_common LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

add_library(pylogger esutil/PyLogger.cpp)
target_link_libraries(pylogger ${BOOST})
install(TARGETS pylogger LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)

if(BUILD_STANDALONE)
  add_executable(pypresso main/pypresso.cpp)
  target_link_libraries(pypresso espresso_common pylogger MPI)
  add_dependencies(pypresso scripts)
  install(TARGETS pypresso RUNTIME DESTINATION bin)
endif(BUILD_STANDALONE)

#we have to build _espresso in a subdir to avoid endless recursion in python
add_subdirectory(pymodule)
