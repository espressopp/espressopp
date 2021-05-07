add_test(NPartSubregion ${Python3_EXECUTABLE} ${PY_COV_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/test_NPartSubregion.py)
set_tests_properties(NPartSubregion PROPERTIES ENVIRONMENT "${TEST_ENV}")
