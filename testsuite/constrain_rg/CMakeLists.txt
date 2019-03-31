add_test(constrain_rg ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_constrain_rg.py)
set_tests_properties(constrain_rg PROPERTIES ENVIRONMENT "${TEST_ENV}")
