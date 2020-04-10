add_test(pure_hdd ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_HeSpaDDAv1.py)
set_tests_properties(pure_hdd PROPERTIES ENVIRONMENT "${TEST_ENV}")
