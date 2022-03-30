if(HPX4ESPP_ENABLED)
    add_test(hpx4espp_HPXRuntime ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_hpx4espp_HPXRuntime.py)
    set_tests_properties(hpx4espp_HPXRuntime PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")

    add_custom_target(hpx4espp_lennard_jones_testdata ALL)
    set(_path ${CMAKE_SOURCE_DIR}/testsuite/vec)
    set(_file lennard_jones_fluid_10000_2048.xyz)
    add_custom_command(TARGET hpx4espp_lennard_jones_testdata COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${_path}/${_file} ${CMAKE_CURRENT_BINARY_DIR}/${_file})
    add_test(hpx4espp_lennard_jones ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_hpx4espp_lennard_jones.py)
    set_tests_properties(hpx4espp_lennard_jones PROPERTIES ENVIRONMENT "${ESP_PY_ENV}")
endif()
