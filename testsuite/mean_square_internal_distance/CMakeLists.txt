add_test(mean_square_internal_distance ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_mean_square_internal_distance.py)
set_tests_properties(mean_square_internal_distance PROPERTIES ENVIRONMENT "${TEST_ENV}")
