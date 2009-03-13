AC_INIT(boost_mpi, 1_35_0)
AC_CONFIG_MACRO_DIR(@macrodir@)
AC_CONFIG_AUX_DIR(build-aux)

AM_INIT_AUTOMAKE([gnu subdir-objects])
AC_LANG([C++])

AXES_IMPORT_EXTRA_PREFIXES

AC_PROG_LIBTOOL
AC_PROG_CXX

AXES_MPI
AM_CONDITIONAL(HAVE_MPI,[test "x$axes_cv_mpi" != "xno"])

AXES_PYTHON
AM_CONDITIONAL(HAVE_PYTHON,[test "x$axes_cv_python" != "xno"])
AM_CONDITIONAL(HAVE_PYTHON_MPI,[test "x$axes_cv_python" != "xno" \
                                     -a "x$axes_cv_mpi" != "xno" ])

pythondir=${libdir}/$axes_python_version
AC_SUBST(pythondir)

dnl mainly to divert the defines of automake - e.g. VERSION collides with
dnl use in boost
AC_CONFIG_HEADERS([acconfig.hpp])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT

if test "x$axes_cv_mpi" = "xno"; then
	AC_MSG_WARN([
***************************************************************
*
* compiler does not seem to be MPI capable, not build boost.mpi
*
***************************************************************
])
fi

if test "x$axes_cv_python" = "xno"; then
	AC_MSG_WARN([
***************************************************************
*
* Python not find, won't build python and python-mpi
*
***************************************************************
])
fi
