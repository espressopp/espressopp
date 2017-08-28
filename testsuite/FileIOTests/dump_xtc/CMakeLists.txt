if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(dump_xtc_testdata ALL)
  foreach(_file expected_with_compression.xtc expected_without_compression.xtc) 
    add_custom_command(TARGET dump_xtc_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(dump_xtc ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dump_xtc.py)
set_tests_properties(dump_xtc PROPERTIES ENVIRONMENT "${TEST_ENV}")
