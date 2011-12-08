######################################
# Create hgversion.hpp file          #
######################################
set(CMAKE_MODULE_PATH ${TOP_SOURCE_DIR}/cmake)
find_package(Mercurial)
if (MERCURIAL_FOUND AND IS_DIRECTORY ${TOP_SOURCE_DIR}/.hg)
  MERCURIAL_HG_INFO(${TOP_SOURCE_DIR} ESPP)
  MESSAGE("Current revision is ${ESPP_HG_ID}")
  set (ESPP_HG_ID ${ESPP_HG_ID})
else (MERCURIAL_FOUND AND IS_DIRECTORY ${TOP_SOURCE_DIR}/.hg)
  set (ESPP_HG_ID)
endif (MERCURIAL_FOUND AND IS_DIRECTORY ${TOP_SOURCE_DIR}/.hg)
set (HG_HEADER "hgversion.hpp")
set (NEW_HG_HEADER "new_hgversion.hpp")
file(WRITE ${NEW_HG_HEADER} "static const std::string hgversion = \"${ESPP_HG_ID}\";\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different ${NEW_HG_HEADER} ${HG_HEADER})
execute_process(COMMAND ${CMAKE_COMMAND} -E remove ${NEW_HG_HEADER})

