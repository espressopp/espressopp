######################################
# System documentation (Doxygen)     #
######################################

find_package(Doxygen QUIET)

if (DOXYGEN_FOUND)
    doxygen_add_docs(doxygen ${PROJECT_SOURCE_DIR}/src)
else (DOXYGEN_FOUND)
    message(WARNING "Not building system documentation because Doxygen not found.")
endif (DOXYGEN_FOUND)

######################################
# User documentation (Sphinx)     #
######################################

find_program(SPHINX_EXECUTABLE
        NAMES sphinx-build
        DOC "Sphinx documentation generation tool (http://www.sphinx-doc.org/)"
        )

if (SPHINX_EXECUTABLE)
    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/sphinx)
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

    add_custom_target(sphinx
            COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH} ${SPHINX_EXECUTABLE} -b html ${SPHINX_SOURCE} ${SPHINX_BUILD}/html
            DEPENDS ${SPHINX_SOURCE}/conf.py _espressopp
            )

    add_custom_target(sphinx-pdf
            COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_BINARY_DIR}:$ENV{PYTHONPATH} ${SPHINX_EXECUTABLE} -b latex ${SPHINX_SOURCE} ${SPHINX_BUILD}/pdf
            COMMAND make -C ${SPHINX_BUILD}/pdf
            DEPENDS ${SPHINX_SOURCE}/conf.py _espressopp
            )
else (SPHINX_EXECUTABLE)
    message(WARNING "Not building user documentation because Sphinx not found.")
endif (SPHINX_EXECUTABLE)
