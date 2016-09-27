add_test(LangevinThermostat ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_LangevinThermostat.py)
set_tests_properties(LangevinThermostat PROPERTIES ENVIRONMENT "${TEST_ENV}")
