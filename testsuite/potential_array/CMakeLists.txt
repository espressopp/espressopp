add_test(potential_array ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/potential_array_test.py)
set_tests_properties(potential_array PROPERTIES ENVIRONMENT "${TEST_ENV}")
