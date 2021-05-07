add_test(minimize_energy ${Python3_EXECUTABLE} ${PY_COV_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/test_minimize_energy.py)
set_tests_properties(minimize_energy PROPERTIES ENVIRONMENT "${TEST_ENV}")
