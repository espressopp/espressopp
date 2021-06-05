add_test(dump_gromacs_adress ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dump_gromacs_adress.py)
set_tests_properties(dump_gromacs_adress PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
