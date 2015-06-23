add_test(h5md_writer ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/h5md_writer.py)
set_tests_properties(h5md_writer PROPERTIES ENVIRONMENT "${TEST_ENV}")
