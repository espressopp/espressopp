add_test(dihedral_harmonic ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_dihedralharmonic.py)
set_tests_properties(dihedral_harmonic PROPERTIES ENVIRONMENT "${TEST_ENV}")
