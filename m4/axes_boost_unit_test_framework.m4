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
#   And sets the ac variables:
#
#     axes_cv_boost_unit_test_framework
#     axes_cv_boost_lib_suffix
#     axes_cv_boost_lib_path
#
# LAST MODIFICATION
#
#   2008-08-15
#
# COPYLEFT
#
#   Copyright (c) 2008 Axel Arnold <axel.arnold@scai.fraunhofer.de>
#
#   based on AX_BOOST_UNIT_TEST_FRAMEWORK, copylefted as follows:
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

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
])
