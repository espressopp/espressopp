AC_INIT([pmi],[0.9])

AC_CONFIG_MACRO_DIR(build-aux/macros)
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_SRCDIR([pmi/pmi.hpp])
AM_INIT_AUTOMAKE()

AXES_IMPORT_EXTRA_PREFIXES

AC_PROG_CXX
AC_LANG([C++])

AC_PROG_LIBTOOL

AXES_MPI
if test "x$axes_cv_mpi" = "xno"; then
   AC_MSG_ERROR([could not compile an MPI program])
fi

AC_MSG_CHECKING([whether to optimize the code])
AC_ARG_WITH([optimize],
   AS_HELP_STRING([--with-optimize], 
     [optimize the code (deactivates logging and control messages, default: no)]),
   [pmi_optimize=$withval], [pmi_optimize="no"])
AC_MSG_RESULT($pmi_optimize)
if test "x$pmi_optimize" = "xyes"; then
   AC_DEFINE(PMI_OPTIMIZE,1,[define if PMI shall be optimized])
   LOG4ESPP_LIBS=""
else
   AXES_LOG4ESPP
fi

CPPFLAGS='-I $(top_srcdir) '"$CPPFLAGS"
LIBS="$LOG4ESPP_LIBS $MPI_LIBS $LIBS"

AC_CONFIG_HEADERS([acconfig.hpp])
AC_CONFIG_FILES([Makefile \
		pmi/Makefile \
		pmi/mpitransmit/Makefile \
		examples/Makefile \
		examples/basic/Makefile \
		examples/proxy/Makefile \
		])

AC_OUTPUT
	
