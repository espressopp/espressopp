add_test(pickle_potential ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pickle_potential.py)
set_tests_properties(pickle_potential PROPERTIES ENVIRONMENT "${TEST_ENV}")
add_test(testwarmup ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testwarmup.py)
set_tests_properties(testwarmup PROPERTIES ENVIRONMENT "${TEST_ENV}")
