add_test(minimize_energy ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_minimize_energy.py)
set_tests_properties(minimize_energy PROPERTIES ENVIRONMENT "${TEST_ENV}")
