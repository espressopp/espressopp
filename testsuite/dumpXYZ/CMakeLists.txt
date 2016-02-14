if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(dumpXYZ_testdata ALL)
  foreach(_file precomputed_xyz.xyz precomputed_comment_line_differ_xyz.xyz test_dumpXYZ_type_not_hardcoded.xyz) 
    add_custom_command(TARGET dumpXYZ_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(dumpXYZ ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dumpXYZ.py)
set_tests_properties(dumpXYZ PROPERTIES ENVIRONMENT "${TEST_ENV}")
