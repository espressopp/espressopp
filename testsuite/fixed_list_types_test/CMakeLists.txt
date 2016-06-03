add_test(fixed_list_types_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fixed_list_types_test.py)
set_tests_properties(fixed_list_types_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
