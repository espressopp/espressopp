add_test(DPDThermostat ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_DPDThermostat.py)
set_tests_properties(DPDThermostat PROPERTIES ENVIRONMENT "${TEST_ENV}")
