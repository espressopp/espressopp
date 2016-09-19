add_test(RDFatomistic ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_RDFatomistic.py)
set_tests_properties(RDFatomistic PROPERTIES ENVIRONMENT "${TEST_ENV}")
