add_test(adress_water ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/system_test.py)
set_tests_properties(adress_water PROPERTIES ENVIRONMENT "${TEST_ENV}")
