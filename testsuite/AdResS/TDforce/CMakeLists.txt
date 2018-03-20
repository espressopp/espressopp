if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(TDforce_testdata ALL)
  foreach(_file table_tf.tab table_tf2.tab)
    add_custom_command(TARGET TDforce_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(TDforce ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_TDforce.py)
set_tests_properties(TDforce PROPERTIES ENVIRONMENT "${TEST_ENV}")
