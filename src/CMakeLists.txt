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

add_custom_target(scripts ALL COMMENT "Copying scripts to builddir: *.py")
foreach(FILENAME ${PYTHON_SCRIPTS})
  set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/${FILENAME}")
  set(DST "${CMAKE_CURRENT_BINARY_DIR}/${FILENAME}")
  if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
    add_custom_command(TARGET scripts COMMAND ${CMAKE_COMMAND} -E copy ${SRC} ${DST})
  endif(NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  install(FILES ${DST} DESTINATION ${PYTHON_INSTDIR}/espresso)
endforeach(FILENAME)

#remove this two lines if autotools are gone
add_custom_command(TARGET scripts COMMAND ${CMAKE_COMMAND} -E touch 
  ${CMAKE_CURRENT_BINARY_DIR}/main/_setupPath.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/main/_setupPath.py DESTINATION ${PYTHON_INSTDIR}/espresso)

#small trick to make PYTHONPATH easy
add_custom_command(TARGET scripts COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_BINARY_DIR}/contrib/mpi4py/MPI.so MPI.so)
add_custom_command(TARGET scripts COMMAND ${CMAKE_COMMAND} -E create_symlink
  ${CMAKE_CURRENT_BINARY_DIR} espresso)

list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES})

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.hpp.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/acconfig.hpp)

set (PYTHON_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (BIN_DIR "${CMAKE_CURRENT_BINARY_DIR}")
set (LIB_DIR "${CMAKE_CURRENT_BINARY_DIR}")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ESPRC.cmakein
  ${CMAKE_BINARY_DIR}/ESPRC @ONLY)
set (PYTHON_DIR "${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR}/espresso")
set (BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin")
set (LIB_DIR "${CMAKE_INSTALL_PREFIX}/lib")
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/ESPRC.cmakein
  ${CMAKE_CURRENT_BINARY_DIR}/ESPRC @ONLY)
install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/ESPRC DESTINATION bin)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

if (EXTERNAL_BOOST)
  set (BOOST_LIB ${Boost_LIBRARIES})
else(EXTERNAL_BOOST)
  set (BOOST espresso_boost)
endif(EXTERNAL_BOOST)

add_library(espresso_common ${ESPRESSO_SOURCES})
target_link_libraries(espresso_common ${BOOST} MPI ${PYTHON_LIBRARIES} ${MPI_LIBRARIES})
install(TARGETS espresso_common LIBRARY DESTINATION lib)

add_executable(pypresso main/pypresso.cpp)
target_link_libraries(pypresso espresso_common)
add_dependencies(pypresso scripts)
install(TARGETS pypresso RUNTIME DESTINATION bin)

add_library(_espresso main/_espresso.cpp)
target_link_libraries(_espresso espresso_common)
add_dependencies(_espresso scripts)
set_target_properties(_espresso PROPERTIES PREFIX "")
install(TARGETS _espresso LIBRARY DESTINATION ${PYTHON_INSTDIR})
