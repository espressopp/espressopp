add_test(isfixed_property ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_isfixed.py)
set_tests_properties(isfixed_property PROPERTIES ENVIRONMENT "${TEST_ENV}")
