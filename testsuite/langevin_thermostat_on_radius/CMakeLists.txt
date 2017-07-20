add_test(langevin_thermostat_on_radius ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_langevin_thermostat_on_radius.py)
set_tests_properties(langevin_thermostat_on_radius PROPERTIES ENVIRONMENT "${TEST_ENV}")
