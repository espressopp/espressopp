if(HPX4ESPP_ENABLED)
    file(GLOB_RECURSE HPX4ESPP_SOURCES *.cpp)
    add_library(_hpx4espp_component ${HPX4ESPP_SOURCES})
    target_compile_definitions(_hpx4espp_component PRIVATE -DHPX4ESPP_ENABLED=1)
    target_link_libraries(_hpx4espp_component PUBLIC HPX::hpx PRIVATE Boost::mpi Boost::python${PYTHON_VERSION_NO_DOT})
else()
    set(HPX4ESPP_SOURCES bindings.cpp utils/HPXUtils.cpp)
    add_library(_hpx4espp_component ${HPX4ESPP_SOURCES})
    target_compile_definitions(_hpx4espp_component PRIVATE -DHPX4ESPP_ENABLED=0)
    target_link_libraries(_hpx4espp_component PRIVATE Boost::python${PYTHON_VERSION_NO_DOT})
endif()

target_link_libraries(_hpx4espp_component PRIVATE Python3::Python)
target_include_directories(_hpx4espp_component PRIVATE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src/include>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/src>)
target_link_libraries(_hpx4espp_component PRIVATE MPI::MPI_CXX)
