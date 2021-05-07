add_test(cap_force ${Python3_EXECUTABLE} ${PY_COV_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/test_cap_force.py)
set_tests_properties(cap_force PROPERTIES ENVIRONMENT "${TEST_ENV}")
