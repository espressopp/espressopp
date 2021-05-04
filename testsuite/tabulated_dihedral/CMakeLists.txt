add_test(tabulated_dihedral ${Python3_EXECUTABLE} ${PY_COV_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/tabulated_dihedral.py)
set_tests_properties(tabulated_dihedral PROPERTIES ENVIRONMENT "${TEST_ENV}")
