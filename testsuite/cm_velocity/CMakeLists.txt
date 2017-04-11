add_test(cm_vel ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_CMVelocity.py)
set_tests_properties(cm_vel PROPERTIES ENVIRONMENT "${TEST_ENV}")
