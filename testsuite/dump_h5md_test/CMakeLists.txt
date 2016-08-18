add_test(dump_h5md ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dump_h5md.py)
set_tests_properties(dump_h5md PROPERTIES ENVIRONMENT "${TEST_ENV}")
