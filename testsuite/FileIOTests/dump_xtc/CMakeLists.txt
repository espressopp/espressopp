add_test(dump_xtc ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dump_xtc.py)
set_tests_properties(dump_xtc PROPERTIES ENVIRONMENT "${TEST_ENV}")
