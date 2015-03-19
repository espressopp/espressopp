if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(polymer_melt_tabulated_testdata ALL)
  foreach(_file polymer_melt.start result.dat) 
    add_custom_command(TARGET polymer_melt_tabulated_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(polymer_melt_tabulated ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polymer_melt_tabulated.py)
set_tests_properties(polymer_melt_tabulated PROPERTIES ENVIRONMENT "${TEST_ENV}")
