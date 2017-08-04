add_test(constrain_rg ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_constrain_rg.py)
set_tests_properties(constrain_rg PROPERTIES ENVIRONMENT "${TEST_ENV}")
