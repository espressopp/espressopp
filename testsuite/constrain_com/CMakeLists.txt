add_test(constrain_com ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_constrain_com.py)
set_tests_properties(constrain_com PROPERTIES ENVIRONMENT "${TEST_ENV}")
