add_test(sub_ens_bond ${Python3_EXECUTABLE} ${PY_COV_OPTS} ${CMAKE_CURRENT_SOURCE_DIR}/sub_ens_bond.py)
set_tests_properties(sub_ens_bond PROPERTIES ENVIRONMENT "${TEST_ENV}")
