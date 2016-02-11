add_test(dumpXYZ ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dumpXYZ.py)
set_tests_properties(dumpXYZ PROPERTIES ENVIRONMENT "${TEST_ENV}")
add_test(difftest "bash" ${CMAKE_CURRENT_SOURCE_DIR}/run_diff.sh precomputed_xyz.xyz precomputed_comment_line_differ_xyz.xyz test_dumpXYZ_type_not_hardcoded.xyz)
set_tests_properties(difftest PROPERTIES ENVIRONMENT "${TEST_ENV}")
