add_test(HAdResS ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_HAdResS.py)
set_tests_properties(HAdResS PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
