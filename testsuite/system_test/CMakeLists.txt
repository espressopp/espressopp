add_test(system_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/system_test.py)
set_tests_properties(system_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
