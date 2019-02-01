if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  add_custom_target(PIAdResS_testdata ALL)
  foreach(_file input_test.xyz FEC_H.dat FEC_O.dat ThdForce_H.dat ThdForce_O.dat tableESP_angle.dat tableESP_bondHH.dat tableESP_bondOH.dat tableESP_HW_HW.dat tableESP_HW_OW.dat tableESP_OW_OW.dat)
    add_custom_command(TARGET PIAdResS_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
  endforeach()
endif()
add_test(PIAdResS ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_PIAdResS.py)
set_tests_properties(PIAdResS PROPERTIES ENVIRONMENT "${TEST_ENV}")
