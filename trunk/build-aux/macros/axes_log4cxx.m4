#
# SYNOPSIS
#
#   AXES_LOG4CXX
#
# DESCRIPTION
#
#   Test for the Apache log4cxx library
#
#   This macro defines:
#
#     HAVE_LOG4CXX
#
# LAST MODIFICATION
#
#   2008-09-10
#
# COPYLEFT
#
#   Copyright (c) 2008 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_LOG4CXX], [
AC_ARG_WITH([log4cxx],
    AS_HELP_STRING([--with-log4cxx@<:@=name@:>@], [use log4cxx (default is yes) - it is possible to specify the exact name of the library (optional)]),
    [
        if test "$withval" = "no"; then
            want_log4cxx="no"
        elif test "$withval" = "yes"; then
            want_log4cxx="yes"
            axes_log4cxx_lib="log4cxx"
        else
            want_log4cxx="yes"
            axes_log4cxx_lib="$withval"
        fi
    ],
    [ want_log4cxx="yes"
      axes_log4cxx_lib="log4cxx"
    ])
if test "x$want_log4cxx" = "xyes"; then
    AC_CACHE_CHECK([for log4cxx library],axes_cv_log4cxx,[
	axes_log4cxx_saved_libs="$LIBS"

	LIBS="-l$axes_log4cxx_lib $axes_log4cxx_saved_libs"
        AC_LANG_PUSH([C++])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[
            @%:@include <log4cxx/logger.h>
            using namespace log4cxx;
        ]], [[
            LoggerPtr rootLogger = Logger::getRootLogger();
            LOG4CXX_DEBUG(rootLogger, "debug message");
        ]])],[ axes_cv_log4cxx="$axes_log4cxx_lib" ],[ axes_cv_log4cxx="no" ])
        AC_LANG_POP([C++])

	LIBS="$axes_log4cxx_saved_libs"
    ])

    if test "x$axes_cv_log4cxx" != "xno"; then
	AC_DEFINE([HAVE_LOG4CXX],[1],[whether log4cxx is available])
    fi
else
	axes_cv_log4cxx=no
fi
])
