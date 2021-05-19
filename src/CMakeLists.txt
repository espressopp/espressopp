file(GLOB_RECURSE ESPRESSO_SOURCES *.cpp)
file(GLOB_RECURSE NOT_ESPRESSO_SOURCES esutil/GenLogger.cpp)
file(GLOB_RECURSE DUMP_XTC_SOURCE io/DumpXTC.cpp io/DumpXTCAdress.cpp)

list(REMOVE_ITEM ESPRESSO_SOURCES ${NOT_ESPRESSO_SOURCES} ${DUMP_XTC_SOURCE})

add_custom_target(gitversion COMMAND ${CMAKE_COMMAND} -DTOP_SOURCE_DIR="${CMAKE_SOURCE_DIR}" -P ${ESPP_CMAKE_DIR}/gitversion.cmake)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/acconfig.hpp.cmakein
        ${CMAKE_CURRENT_BINARY_DIR}/acconfig.hpp)

#build .env for pipenv
#CMAKE_BINARY_DIR for espressopp, CMAKE_CURRENT_BINARY_DIR for for _espressopps.so and contrib for MPI.so
if (EXTERNAL_MPI4PY)
    set (PYTHON_DIR "${CMAKE_BINARY_DIR}")
else (EXTERNAL_MPI4PY)
    set (PYTHON_DIR "${CMAKE_BINARY_DIR}:${CMAKE_BINARY_DIR}/contrib")
endif (EXTERNAL_MPI4PY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/pipenv-env.cmakein
        ${CMAKE_BINARY_DIR}/.env @ONLY)

add_library(_espressopp ${ESPRESSO_SOURCES})

if(ESPP_WERROR)
    target_compile_options(_espressopp PRIVATE "-Werror")
endif()
# compiler specific configurations
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    if(ESPP_WALL)
        target_compile_options(_espressopp PRIVATE -Wall)
    endif()
    if(ESPP_LOCAL_ARCHITECTURE)
        target_compile_options(_espressopp PRIVATE "-march=native")
    endif()
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    if(ESPP_LOCAL_ARCHITECTURE)
        target_compile_options(_espressopp PRIVATE "-march=native")
    endif()
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    target_compile_options(_espressopp PRIVATE "-restrict")
    if(ESPP_LOCAL_ARCHITECTURE)
        target_compile_options(_espressopp PRIVATE "-xHost")
    endif()
    if(ESPP_VEC_REPORT)
        target_compile_options(_espressopp PRIVATE "-qopt-report-phase=vec,loop -qopt-report=5")
    endif()
endif()
get_target_property(ESPP_CFLAGS _espressopp COMPILE_OPTIONS)
message(STATUS "_espressopp compiler flags are: ${ESPP_CFLAGS}")

target_include_directories(_espressopp PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src> $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/include>
     $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>)

target_link_libraries(_espressopp PUBLIC Boost::mpi Boost::serialization Boost::system Boost::filesystem Boost::python${PYTHON_VERSION_NO_DOT} Boost::numpy${PYTHON_VERSION_NO_DOT})
target_link_libraries(_espressopp PUBLIC Python3::Python)
target_link_libraries(_espressopp PUBLIC MPI::MPI_CXX)
target_link_libraries(_espressopp PRIVATE FFTW3::fftw3)

if(WITH_XTC)
    target_add_definitions(_espressopp PRIVATE -DHAS_GROMACS)
    target_link_libraries(_espressopp PRIVATE Gromacs::libgromacs)
    target_sources(_espressopp PRIVATE ${DUMP_XTC_SOURCE})
endif()
#python libs have not prefix (default would be 'lib')
set_target_properties(_espressopp PROPERTIES PREFIX "" SUFFIX ".so" LIBRARY_OUTPUT_DIRECTORY ..)
add_dependencies(_espressopp gitversion)

install(TARGETS _espressopp LIBRARY DESTINATION ${PYTHON_INSTDIR} ARCHIVE DESTINATION ${PYTHON_INSTDIR})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ DESTINATION ${PYTHON_INSTDIR}/espressopp
        FILES_MATCHING PATTERN "*.py"
        PATTERN "CMakeFiles" EXCLUDE)

