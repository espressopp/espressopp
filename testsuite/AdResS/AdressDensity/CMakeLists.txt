add_test(AdressDensity ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_AdressDensity.py)
set_tests_properties(AdressDensity PROPERTIES ENVIRONMENT "${TEST_ENV}")
