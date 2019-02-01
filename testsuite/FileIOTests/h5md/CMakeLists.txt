add_test(dump_h5md_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dump_h5md_test.py)
add_test(dump_topology_test ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/dump_topology_test.py)
set_tests_properties(dump_h5md_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
set_tests_properties(dump_topology_test PROPERTIES ENVIRONMENT "${TEST_ENV}")
