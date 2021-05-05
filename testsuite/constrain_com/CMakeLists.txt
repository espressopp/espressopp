add_test(constrain_com ${Python3_EXECUTABLE} ${PY_COV_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/test_constrain_com.py)
set_tests_properties(constrain_com PROPERTIES ENVIRONMENT "${TEST_ENV}")
