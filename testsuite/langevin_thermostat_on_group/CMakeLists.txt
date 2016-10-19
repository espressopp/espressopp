add_test(langevin_thermostat_on_group ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/langevin_thermostat_on_group.py)
set_tests_properties(langevin_thermostat_on_group PROPERTIES ENVIRONMENT "${TEST_ENV}")
