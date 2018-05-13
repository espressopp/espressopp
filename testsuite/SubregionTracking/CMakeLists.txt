add_test(SubregionTracking ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_SubregionTracking.py)
set_tests_properties(SubregionTracking PROPERTIES ENVIRONMENT "${TEST_ENV}")
