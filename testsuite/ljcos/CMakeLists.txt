add_test(ljcos ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_ljcos.py)
set_tests_properties(ljcos PROPERTIES ENVIRONMENT "${TEST_ENV}")
