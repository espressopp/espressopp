add_test(sub_ens_bond ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/sub_ens_bond.py)
set_tests_properties(sub_ens_bond PROPERTIES ENVIRONMENT "${TEST_ENV}")
