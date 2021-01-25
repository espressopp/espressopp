add_test(pure_hdd ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_HeSpaDDAv1p1.py)
set_tests_properties(pure_hdd PROPERTIES ENVIRONMENT "${TEST_ENV}")
