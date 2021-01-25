add_test(dump_gromacs ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dump_gromacs.py)
set_tests_properties(dump_gromacs PROPERTIES ENVIRONMENT "${TEST_ENV}")
