if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(FreeEnergyCompensation_testdata ALL)
  foreach(_file table_fec.tab)
    add_custom_command(TARGET FreeEnergyCompensation_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(FreeEnergyCompensation ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_FreeEnergyCompensation.py)
set_tests_properties(FreeEnergyCompensation PROPERTIES ENVIRONMENT "${TEST_ENV}")
