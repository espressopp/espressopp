add_test(thermo_integration ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_thermo_integration.py)
set_tests_properties(thermo_integration PROPERTIES ENVIRONMENT "${TEST_ENV}")
