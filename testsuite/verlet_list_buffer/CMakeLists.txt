if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(verlet_list_buffer_testdata ALL)
  set(_file lennard_jones_fluid_10000.xyz)
  add_custom_command(TARGET verlet_list_buffer_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endif()
add_test(verlet_list_buffer ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verlet_list_buffer.py)
set_tests_properties(verlet_list_buffer PROPERTIES ENVIRONMENT "${TEST_ENV}")
