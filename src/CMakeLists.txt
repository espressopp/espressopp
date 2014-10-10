file(GLOB_RECURSE ESPRESSO_SOURCES *.cpp)
file(GLOB_RECURSE NOT_ESPRESSO_SOURCES
  main/_espresso.cpp
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

add_custom_target(hgversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/hgversion.cmake)

list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.hpp.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/acconfig.hpp)

#build ESPRC
#CMAKE_BINARY_DIR for espresso, CMAKE_CURRENT_BINARY_DIR for for _espressos.so and contrib for MPI.so
if (EXTERNAL_MPI4PY)
  set (PYTHON_DIR "${CMAKE_BINARY_DIR}:${CMAKE_CURRENT_BINARY_DIR}/pymodule")
else (EXTERNAL_MPI4PY)
  set (PYTHON_DIR "${CMAKE_BINARY_DIR}:${CMAKE_CURRENT_BINARY_DIR}/pymodule:${CMAKE_BINARY_DIR}/contrib")
endif (EXTERNAL_MPI4PY)
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
option(WITH_RC_FILES "Install Espresso rc file, no need when installing under /usr" ON)
if (WITH_RC_FILES AND BUILD_PYTHON_MODULE)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ESPRC DESTINATION bin)
endif (WITH_RC_FILES AND BUILD_PYTHON_MODULE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if (EXTERNAL_BOOST)
  set (BOOST ${Boost_LIBRARIES} )
else(EXTERNAL_BOOST)
  set (BOOST espresso_boost)
endif(EXTERNAL_BOOST)

add_library(espresso_common ${ESPRESSO_SOURCES})
target_link_libraries(espresso_common pylogger ${BOOST} ${PYTHON_LIBRARIES}
                                      ${VAMPIRTRACE_LIBRARIES} ${MPI_LIBRARIES} ${FFTW3_LIBRARIES})
install(TARGETS espresso_common LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION lib)

add_library(pylogger esutil/PyLogger.cpp)
target_link_libraries(pylogger ${BOOST})
install(TARGETS pylogger LIBRARY DESTINATION ${LIB} ARCHIVE DESTINATION lib)

#we have to build _espresso in a subdir to avoid endless recursion in python
add_subdirectory(pymodule)


add_dependencies(espresso_common hgversion)
