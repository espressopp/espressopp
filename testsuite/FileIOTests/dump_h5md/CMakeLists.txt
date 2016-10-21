add_test(dump_h5md_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dump_h5md_test.py)
set_tests_properties(dump_h5md_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
