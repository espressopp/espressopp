add_test(io ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_python_tools_io.py)
set_tests_properties(io PROPERTIES ENVIRONMENT "${TEST_ENV}")
