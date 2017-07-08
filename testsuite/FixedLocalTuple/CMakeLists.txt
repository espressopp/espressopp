add_test(fixedlocaltuple ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/TestFixedLocalTupleList.py)
set_tests_properties(fixedlocaltuple PROPERTIES ENVIRONMENT "${TEST_ENV}")
