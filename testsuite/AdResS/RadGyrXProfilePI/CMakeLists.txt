add_test(RadGyrXProfilePI ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_RadGyrXProfilePI.py)
set_tests_properties(RadGyrXProfilePI PROPERTIES ENVIRONMENT "${TEST_ENV}")
