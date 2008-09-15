#
# SYNOPSIS
#
#   AXES_BOOST_SERIALIZATION
#
# DESCRIPTION
#
#   Test for the serialization library from the Boost C++ libraries.
#   The macro requires a preceding call to AXES_BOOST_BASE.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_SERIALIZATION_LIBS)
#
#   And defines:
#
#     HAVE_BOOST_SERIALIZATION
#
#   And sets the ac variables:
#
#     axes_cv_boost_serialization
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
#   based on AX_BOOST_SERIALIZATION, copylefted as follows:
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_BOOST_SERIALIZATION],
[
AC_REQUIRE([AXES_BOOST_BASE])
AXES_BOOST_LIB(serialization, SERIALIZATION, AC_LANG_PROGRAM([[
        @%:@include <boost/archive/basic_archive.hpp>
]], [[
        boost::archive::ARCHIVE_VERSION();
        return 0;
]]))
])
