#
# SYNOPSIS
#
#   AXES_BOOST_UNIT_TEST_FRAMEWORK
#
# DESCRIPTION
#
#   Test for the unit test framework from the Boost C++ libraries.
#   The macro requires a preceding call to AXES_BOOST_BASE.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIBS)
#
#   And defines:
#
#     HAVE_BOOST_UNIT_TEST_FRAMEWORK
#
#   If the dynamic unit test framework library is used, it also defines
#
#     BOOST_TEST_DYN_LINK
#
#   And sets the ac variables:
#
#     axes_cv_boost_unit_test_framework
#     axes_cv_boost_lib_suffix
#     axes_cv_boost_lib_path
#
# LAST MODIFICATION
#
#   2009-02-20
#
# COPYLEFT
#
#   Copyright (c) 2008 Axel Arnold <axel.arnold@scai.fraunhofer.de>
#   Copyright (c) 2009 Olaf Lenz <lenzo@mpip-mainz.mpg.de>
#
#   based on AX_BOOST_UNIT_TEST_FRAMEWORK, copylefted as follows:
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

# first test whether the library is there and can be linked
AC_DEFUN([AXES_BOOST_UNIT_TEST_FRAMEWORK],
[
AC_REQUIRE([AXES_BOOST_BASE])
AXES_BOOST_LIB(unit_test_framework, UNIT_TEST_FRAMEWORK, AC_LANG_PROGRAM([[
        @%:@include <boost/test/unit_test.hpp>
]], [[

        using boost::unit_test_framework::test_suite;
        using boost::unit_test_framework::test_case;

        test_suite* test= BOOST_TEST_SUITE( "operty Test" );

]]))
	
# now link a test program and test whether we have the static or dynamic version

# save current flags
axes_boost_unit_test_framework_cppflags_saved="$CPPFLAGS"
axes_boost_unit_test_framework_ldflags_saved="$LDFLAGS"
axes_boost_unit_test_framework_libs_saved="$LIBS"

# set up the boost flags
CPPFLAGS="$BOOST_CPPFLAGS $CPPFLAGS"
LDFLAGS="$BOOST_LDFLAGS $LDFLAGS"
LIBS="$BOOST_UNIT_TEST_FRAMEWORK_LIBS $LIBS"

# test link a program using the dynamic linking
AC_MSG_CHECKING([whether a test can be dynamically linked])
AC_LANG_PUSH([C++])
AC_LINK_IFELSE(AC_LANG_SOURCE([[
        @%:@define BOOST_TEST_DYN_LINK
        @%:@define BOOST_TEST_MODULE conftest
        @%:@include <boost/test/unit_test.hpp>
        BOOST_AUTO_TEST_CASE(conftestcase) {
          BOOST_CHECK(true);
        }
]]), 
axes_boost_unit_test_framework_succeeded=yes,
axes_boost_unit_test_framework_succeeded=no)
AC_LANG_POP([C++])
AC_MSG_RESULT([$axes_boost_unit_test_framework_succeeded])

# define BOOST_TEST_DYN_LINK if required
if test "x$axes_boost_unit_test_framework_succeeded" = "xyes"; then
  AC_DEFINE([BOOST_TEST_DYN_LINK], 1, 
  [whether to link the unit test framework library dynamically])
fi

# restore previous flags
CPPFLAGS="$axes_boost_unit_test_framework_cppflags_saved"
LDFLAGS="$axes_boost_unit_test_framework_ldflags_saved"
LIBS="$axes_boost_unit_test_framework_libs_saved"

])
