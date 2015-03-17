file(GLOB_RECURSE ESPRESSO_SOURCES *.cpp)
file(GLOB_RECURSE NOT_ESPRESSO_SOURCES
  esutil/GenLogger.cpp
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

file(GLOB_RECURSE PYTHON_SCRIPTS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)

add_custom_target(gitversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${CMAKE_MODULE_PATH}/gitversion.cmake)

list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.hpp.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/acconfig.hpp)

#build ESPRC
#CMAKE_BINARY_DIR for espressopp, CMAKE_CURRENT_BINARY_DIR for for _espressopps.so and contrib for MPI.so
if (EXTERNAL_MPI4PY)
  set (PYTHON_DIR "${CMAKE_BINARY_DIR}")
else (EXTERNAL_MPI4PY)
  set (PYTHON_DIR "${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/contrib")
endif (EXTERNAL_MPI4PY)
#CMAKE_CURRENT_BINARY_DIR for libespressopp_common.so and contib for libespressopp_boost
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ESPRC.cmakein
  ${CMAKE_BINARY_DIR}/ESPRC @ONLY)

#after install ESPRC
set (PYTHON_DIR "${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ESPRC.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/ESPRC @ONLY)
option(WITH_RC_FILES "Install Espresso rc file, no need when installing under /usr" ON)
if (WITH_RC_FILES)
  install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ESPRC DESTINATION bin)
endif (WITH_RC_FILES)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if (EXTERNAL_BOOST)
  set (BOOST ${Boost_LIBRARIES} )
else(EXTERNAL_BOOST)
  set (BOOST espressopp_boost)
endif(EXTERNAL_BOOST)

add_library(_espressopp ${ESPRESSO_SOURCES})
target_link_libraries(_espressopp ${BOOST} ${PYTHON_LIBRARIES} ${MPI_LIBRARIES} ${FFTW3_LIBRARIES} ${VAMPIRTRACE_LIBRARIES})
#python libs have not prefix (default would be 'lib')
set_target_properties(_espressopp PROPERTIES PREFIX "" SUFFIX ".so" LIBRARY_OUTPUT_DIRECTORY ..)
install(TARGETS _espressopp LIBRARY DESTINATION ${PYTHON_INSTDIR})

add_dependencies(_espressopp gitversion)
