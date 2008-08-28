AC_INIT(boost_mpi, 1_35_0)
AC_CONFIG_MACRO_DIR(../../../config)
AC_CONFIG_AUX_DIR(build-aux)

AM_INIT_AUTOMAKE(1.9)
AC_LANG([C++])

AXES_IMPORT_PREFIXES

AC_PROG_LIBTOOL
AC_PROG_CXX

dnl check MPI
AXES_MPI
if test "x$axes_cv_mpi" = "xno"; then
	AC_MSG_ERROR([compiler does not seem to be MPI capable])
fi

dnl look for Boost
AXES_BOOST_BASE([1.35])
AXES_BOOST_SERIALIZATION
if test "x$axes_cv_boost_serialization" = "xno"; then
	AC_MSG_ERROR([could not find the boost serialization library])
fi

if test "$axes_cv_boost_lib_suffix" = ".none."; then
	BOOST_SUFFIX=
else
	BOOST_SUFFIX=$axes_cv_boost_lib_suffix
fi
AC_SUBST(BOOST_SUFFIX)

dnl load our headers for boost.mpi _before_ the system-installed
dnl it can happen that there are the headers in the system, but
dnl not the library. Then, we need to compile against our own headers,
dnl since they can differ (e.g., in version)
CPPFLAGS=" -I\$(top_builddir) $BOOST_CPPFLAGS"
LIBS="\$(BOOST_SERIALIZATION_LIBS)"

dnl check if the boost headers were in a versioned dir, then version ourselves
if test "x$axes_cv_boost" != "x"; then
    p=${axes_cv_boost%/}
    p=${p##*/}
    includedir="$includedir/$p"
fi

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
