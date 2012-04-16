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
  "libs/detail/utf8_codecvt_facet.cpp"
)

list(REMOVE_ITEM BOOST_SOURCES ${NOT_BOOST_SOURCES})

add_library(espresso_boost ${BOOST_SOURCES})
target_link_libraries(espresso_boost ${PYTHON_LIBRARIES} ${VAMPIRTRACE_LIBRARIES} ${MPI_LIBRARIES})
install(TARGETS espresso_boost LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
