add_test(verlet_list_dynamic_exclude_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verlet_list_dynamic_exclude_test.py)
set_tests_properties(verlet_list_dynamic_exclude_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
