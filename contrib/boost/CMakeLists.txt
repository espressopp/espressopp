file(GLOB_RECURSE BOOST_SOURCES libs/*.cpp)
file(GLOB_RECURSE NOT_BOOST_TESTS "libs/test/*")
file(GLOB_RECURSE NOT_BOOST_SOURCES
  "libs/test/src/cpp_main.cpp"
  "libs/test/src/test_main.cpp"
  "libs/serialization/src/basic_text_wiprimitive.cpp"
  "libs/serialization/src/basic_text_woprimitive.cpp"
  "libs/serialization/src/text_wiarchive.cpp"
  "libs/serialization/src/text_woarchive.cpp"
  "libs/serialization/src/utf8_codecvt_facet.cpp"
  "libs/serialization/src/xml_wgrammar.cpp"
  "libs/serialization/src/xml_wiarchive.cpp"
  "libs/serialization/src/xml_woarchive.cpp"
  "libs/serialization/src/codecvt_null.cpp"
  "libs/serialization/src/binary_wiarchive.cpp"
  "libs/serialization/src/binary_woarchive.cpp"
  "libs/serialization/src/shared_ptr_helper.cpp"
  "libs/detail/utf8_codecvt_facet.cpp"
  "libs/atomic/src/lockpool.cpp"
  "libs/thread/src/pthread/once_atomic.cpp"
  "libs/thread/src/pthread/once.cpp"
  "libs/thread/src/pthread/thread.cpp"
)

list(REMOVE_ITEM BOOST_SOURCES ${NOT_BOOST_SOURCES})
set(BOOST_SOURCES_TEST ${BOOST_SOURCES})
# Remove from the main boost library also tests framework.
list(REMOVE_ITEM BOOST_SOURCES ${NOT_BOOST_TESTS})

# Apply patch for Intel Compiler; http://stackoverflow.com/a/41188485
if (NOT EXTERNAL_BOOST AND CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(STATUS "Applying patch on internal Boost because of Intel compiler")
    execute_process(COMMAND patch
        ${CMAKE_CURRENT_LIST_DIR}/boost/serialization/array.hpp
        ${CMAKE_CURRENT_LIST_DIR}/serialization_array_hpp.patch
        RESULT_VARIABLE boost_patch OUTPUT_VARIABLE boost_patch_stdout)
    #message(${boost_patch_stdout})
    if ("${boost_patch}" EQUAL "1" AND NOT "${boost_patch_stdout}" MATCHES "Skipping patch")
        message(FATAL_ERROR ${boost_patch_stdout})
    endif()
endif()

# Static boost library for the main _espressopp.so module.
add_library(espressopp_boost STATIC ${BOOST_SOURCES})
# Shared libboost library for unit tests.
add_library(espressopp_boost_test SHARED ${BOOST_SOURCES_TEST})
target_link_libraries(espressopp_boost ${PYTHON_LIBRARIES} ${VAMPIRTRACE_LIBRARIES} ${MPI_LIBRARIES} ${RT_LIBRARIES})
if( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  set_target_properties(espressopp_boost PROPERTIES COMPILE_FLAGS "-fPIC")
  set_target_properties(espressopp_boost_test PROPERTIES COMPILE_FLAGS "-fPIC")
endif( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
