add_test(python_tools ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_python_tools.py)
set_tests_properties(python_tools PROPERTIES ENVIRONMENT "${TEST_ENV}")
