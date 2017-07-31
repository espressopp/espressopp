add_test(complex ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_python_tools_complex.py)
set_tests_properties(complex PROPERTIES ENVIRONMENT "${TEST_ENV}")
