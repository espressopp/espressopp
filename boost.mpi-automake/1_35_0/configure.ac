AC_INIT(boost_mpi, 1_35_0)
AM_INIT_AUTOMAKE(1.9)
AC_LANG([C++])

dnl look for MPI
ACX_MPI
CXX="$MPICXX"

dnl look for Boost
AX_BOOST_BASE([1.35])
AX_BOOST_SERIALIZATION
if test "x$ax_cv_boost_serialization" != "xyes"; then
	AC_MSG_ERROR([could not find the boost serialization library])
fi

dnl load our headers for boost.mpi _before_ the system-installed
dnl it can happen that there are the headers in the system, but
dnl not the library. Then, we need to compile against our own headers,
dnl since they can differ (e.g., in version)
BOOST_CPPFLAGS=" -I\$(top_builddir) $BOOST_CPPFLAGS"
AC_SUBST(BOOST_CPPFLAGS)

AC_PROG_LIBTOOL
AC_SUBST(LIBTOOL_DEPS)

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile src/Makefile])

AC_OUTPUT
