add_test(RadGyrXProfilePI ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_RadGyrXProfilePI.py)
set_tests_properties(RadGyrXProfilePI PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
