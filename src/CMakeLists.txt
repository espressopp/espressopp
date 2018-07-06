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

if(NOT WITH_XTC)
  file(GLOB_RECURSE DUMP_XTC_SOURCE io/DumpXTC.cpp io/DumpXTCAdress.cpp)
  list(REMOVE_ITEM ESPRESSO_SOURCES ${DUMP_XTC_SOURCE})
endif()

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

add_library(_espressopp ${ESPRESSO_SOURCES})
target_link_libraries(_espressopp ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} ${MPI_LIBRARIES} ${FFTW3_LIBRARIES} ${VAMPIRTRACE_LIBRARIES})
if(WITH_XTC)
  target_link_libraries(_espressopp ${GROMACS_LIBRARIES})
endif()
#python libs have not prefix (default would be 'lib')
set_target_properties(_espressopp PROPERTIES PREFIX "" SUFFIX ".so" LIBRARY_OUTPUT_DIRECTORY ..)
install(TARGETS _espressopp LIBRARY DESTINATION ${PYTHON_INSTDIR} ARCHIVE DESTINATION ${PYTHON_INSTDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION ${PYTHON_INSTDIR}/espressopp
	FILES_MATCHING PATTERN "*.py"
    PATTERN "CMakeFiles" EXCLUDE)

add_dependencies(_espressopp gitversion)

find_program(CLANG_FORMAT "clang-format")
if (CLANG_FORMAT)
    add_custom_target(
            format
            COMMAND ${CLANG_FORMAT}
            -style=file
            -i
            ${ESPRESSO_SOURCES}
    )
endif()
