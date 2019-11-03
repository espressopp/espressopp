if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(potential_array_testdata ALL)
  set(_file lennard_jones_fluid_10000_2048.xyz)
  add_custom_command(TARGET potential_array_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endif()
add_test(potential_array ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/potential_array_test.py)
set_tests_properties(potential_array PROPERTIES ENVIRONMENT "${TEST_ENV}")
