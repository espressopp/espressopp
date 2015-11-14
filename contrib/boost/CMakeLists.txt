file(GLOB_RECURSE BOOST_SOURCES libs/*.cpp)
file(GLOB_RECURSE NOT_BOOST_SOURCES
  libs/test/*
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

add_library(espressopp_boost STATIC ${BOOST_SOURCES})
target_link_libraries(espressopp_boost ${PYTHON_LIBRARIES} ${VAMPIRTRACE_LIBRARIES} ${MPI_LIBRARIES} ${RT_LIBRARIES})
if( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
  set_target_properties(espressopp_boost PROPERTIES COMPILE_FLAGS "-fPIC")
endif( CMAKE_SYSTEM_PROCESSOR STREQUAL "x86_64" )
