add_test(fixed_pair_list_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fixed_pair_list_test.py)
set_tests_properties(fixed_pair_list_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
