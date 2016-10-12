add_test(dump_gromacs_adress ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dump_gromacs_adress.py)
set_tests_properties(dump_gromacs_adress PROPERTIES ENVIRONMENT "${TEST_ENV}")
