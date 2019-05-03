add_test(ForceAdResS ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_AdResS.py)
set_tests_properties(ForceAdResS PROPERTIES ENVIRONMENT "${TEST_ENV}")
