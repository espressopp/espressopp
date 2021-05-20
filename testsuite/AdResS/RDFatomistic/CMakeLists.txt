add_test(RDFatomistic ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_RDFatomistic.py)
set_tests_properties(RDFatomistic PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
