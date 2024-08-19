# Lennard-Jones system
if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(vec_lennard_jones_testdata ALL)
  set(_file lennard_jones_fluid_10000_2048.xyz)
  add_custom_command(TARGET vec_lennard_jones_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endif()
add_test(vec_lennard_jones ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_vec_lennard_jones.py)
set_tests_properties(vec_lennard_jones PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")

# Polymer melt system
if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(vec_polymer_melt_testdata ALL)
  set(_file polymer_melt.lammps)
  add_custom_command(TARGET vec_polymer_melt_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
    ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
endif()
add_test(vec_polymer_melt ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_vec_polymer_melt.py)
set_tests_properties(vec_polymer_melt PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")

# Langevin Thermostat
add_test(vec_langevin_thermostat ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_vec_langevin_thermostat.py)
set_tests_properties(vec_langevin_thermostat PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
