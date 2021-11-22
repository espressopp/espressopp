configure_file(${CMAKE_CURRENT_SOURCE_DIR}/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/. COPYONLY)
add_test(h5md_parallel ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_h5md_parallel.py)
set_tests_properties(h5md_parallel PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
