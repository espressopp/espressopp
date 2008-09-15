#
# SYNOPSIS
#
#   AXES_LOG4CPP
#
# DESCRIPTION
#
#   Test for the log4cpp library
#
#   This macro defines:
#
#     HAVE_LOG4CPP
#
# LAST MODIFICATION
#
#   2008-09-09
#
# COPYLEFT
#
#   Copyright (c) 2008 Thomas Brandes <thomas.brandes@scai.fraunhofer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_LOG4CPP], [
AC_ARG_WITH([log4cpp],
    AS_HELP_STRING([--with-log4cpp@<:@=name@:>@], [use log4cpp (default is yes) - it is possible to specify the exact name of the library (optional)]),
    [
        if test "$withval" = "no"; then
            want_log4cpp="no"
        elif test "$withval" = "yes"; then
            want_log4cpp="yes"
            axes_log4cpp_lib="log4cpp"
        else
            want_log4cpp="yes"
            axes_log4cpp_lib="$withval"
        fi
    ],
    [ want_log4cpp="yes"
      axes_log4cpp_lib="log4cpp"
    ])

if test "x$want_log4cpp" = "xyes"; then
    AC_CACHE_CHECK([for log4cpp library],axes_cv_log4cpp,[
	axes_log4cpp_saved_libs="$LIBS"

	LIBS="-l$axes_log4cpp_lib $axes_log4cpp_saved_libs"

        AC_LANG_PUSH([C++])
        AC_LINK_IFELSE([AC_LANG_PROGRAM([[
            @%:@include <log4cpp/Category.hh>
        ]], [[
            log4cpp::Category& myLogger = log4cpp::Category::getRoot();
            LOG4CPP_DEBUG(myLogger, "debug message");
        ]])],[ axes_cv_log4cpp="$axes_log4cpp_lib" ],[ axes_cv_log4cpp="no" ])
	AC_LANG_POP([C++])

	LIBS="$axes_log4cpp_saved_libs"
    ])

    if test "x$axes_cv_log4cpp" != "xno"; then
	AC_DEFINE([HAVE_LOG4CPP],[1],[whether log4cpp is available])
    fi
else
    axes_cv_log4cpp=no
fi
])
