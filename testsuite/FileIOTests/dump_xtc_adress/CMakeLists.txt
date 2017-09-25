if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(dump_xtc_adress_testdata ALL)
  add_custom_command(TARGET dump_xtc_adress_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/expected.xtc ${CMAKE_CURRENT_BINARY_DIR}/expected.xtc)
endif()
add_test(dump_xtc_adress ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dump_xtc_adress.py)
set_tests_properties(dump_xtc_adress PROPERTIES ENVIRONMENT "${TEST_ENV}")
