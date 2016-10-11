add_test(pure_lb ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_LatticeBoltzmann.py)
set_tests_properties(pure_lb PROPERTIES ENVIRONMENT "${TEST_ENV}")
