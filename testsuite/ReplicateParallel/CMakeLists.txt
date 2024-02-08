if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(ReplicateParallel_testdata ALL)
  set(_file polymer_melt.lammps)
  add_custom_command(TARGET ReplicateParallel_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endif()
message(STATUS "Python3_EXECUTABLE=${Python3_EXECUTABLE}")
add_test(ReplicateParallel ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/polymer_melt.py)
set_tests_properties(ReplicateParallel PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
