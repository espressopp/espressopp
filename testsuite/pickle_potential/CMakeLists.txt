add_test(pickle_potential ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/pickle_potential.py)
set_tests_properties(pickle_potential PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")
add_test(testwarmup ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/testwarmup.py)
set_tests_properties(testwarmup PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH}")
