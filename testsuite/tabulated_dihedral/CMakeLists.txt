add_test(tabulated_dihedral ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tabulated_dihedral.py)
set_tests_properties(tabulated_dihedral PROPERTIES ENVIRONMENT "${TEST_ENV}")
