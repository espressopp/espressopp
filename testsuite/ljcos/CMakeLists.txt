add_test(ljcos ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_ljcos.py)
set_tests_properties(ljcos PROPERTIES ENVIRONMENT "${TEST_ENV}")
