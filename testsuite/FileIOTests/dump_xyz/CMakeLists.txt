add_test(dump_xyz ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dump_xyz.py)
set_tests_properties(dump_xyz PROPERTIES ENVIRONMENT "${TEST_ENV}")
