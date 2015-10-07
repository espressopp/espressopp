if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
 add_custom_target(ewald_eppDeserno_comparison_testdata ALL)
  foreach(_file deserno_ewald.dat ini_struct_deserno.dat)
    add_custom_command(TARGET ewald_eppDeserno_comparison_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(ewald_eppDeserno_comparison ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/ewald_eppDeserno_comparison.py)
set_tests_properties(ewald_eppDeserno_comparison PROPERTIES ENVIRONMENT "${TEST_ENV}")
