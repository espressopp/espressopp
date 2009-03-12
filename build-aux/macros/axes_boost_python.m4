#
# SYNOPSIS
#
#   AXES_BOOST_PYTHON
#
# DESCRIPTION
#
#   Test for the Python bindings library from the Boost C++ libraries.
#   The macro requires a preceding call to AXES_BOOST_BASE and to
#   AX_PYTHON.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_PYTHON_LIBS)
#
#   And defines:
#
#     HAVE_BOOST_PYTHON
#
#   And sets the ac variables:
#
#     axes_cv_boost_python
#     axes_cv_boost_lib_suffix
#     axes_cv_boost_lib_path
#
#   Furthermore, it appends the contents of the variable
#   PYTHON_CPPFLAGS to the variable BOOST_CPPFLAGS.
#
# LAST MODIFICATION
#
#   2009-03-11
#
# COPYLEFT
#
#   Copyright (c) 2009 Olaf Lenz <lenzo@mpip-mainz.mpg.de>
#   Copyright (c) 2008 Axel Arnold <axel.arnold@scai.fraunhofer.de>
#
#   based on AX_BOOST_SERIALIZATION, copylefted as follows:
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_BOOST_PYTHON],
[
AC_REQUIRE([AXES_PYTHON])

if test "x$axes_cv_python_lib" != "xno" && test "x$axes_cv_python_include" != "xno"; then
	dnl temporarily add python flags
	axes_boost_python_saved_cppflags="$CPPFLAGS"
	axes_boost_python_saved_libs="$LIBS"
	CPPFLAGS="$axes_boost_python_saved_cppflags $PYTHON_CPPFLAGS"
	LIBS="$axes_boost_python_saved_libs $PYTHON_LIBS"

	AXES_BOOST_LIB(python, PYTHON, AC_LANG_PROGRAM([[
		 #include <boost/python.hpp>
		 using namespace boost::python;
	]], [[
		boost::python::detail::init_module("test", 0);
	        return 0;
	]]))

	BOOST_CPPFLAGS="$BOOST_CPPFLAGS $PYTHON_CPPFLAGS"
	AC_SUBST(BOOST_CPPFLAGS)
	BOOST_PYTHON_LIBS="$BOOST_PYTHON_LIBS $PYTHON_LIBS"
	AC_SUBST(BOOST_PYTHON_LIBS)
	dnl restore flags
	CPPFLAGS="$axes_boost_python_saved_cppflags"
	LIBS="$axes_boost_python_saved_libs"
else
	axes_cv_boost_lib_python=no
fi
])
