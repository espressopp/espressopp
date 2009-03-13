# DESCRIPTION
#   
#   Chooses the logger library for LOG4ESPP.
#
#   This macro calls:
#
#     AC_SUBST(LOG4ESPP_LIBS)
#     AC_SUBST(LOG4ESPP_USE)
#
#   And sets the ac variables:
#
#     axes_log4espp
#
# LAST MODIFICATION
#
#   2009-03-12
#
# COPYLEFT
#
#   Copyright (c) 2009 Olaf Lenz <lenzo@mpip-mainz.mpg.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_LOG4ESPP],
[
AC_REQUIRE([AXES_LOG4CXX])
AC_REQUIRE([AXES_LOG4CPP])

AC_ARG_WITH([log4espp],
   AS_HELP_STRING([--with-log4espp=LOGLIB], 
[use the specified logger library (default: auto). 
Possible values are "log4cxx", "log4cpp", "generic", "no" and "auto".])
  ,, [with_log4espp="auto"])

AC_MSG_CHECKING([which logger library to use])

dnl Test which logger to take
if test "x$with_log4espp" = "xauto" || test "x$with_log4espp" = "xyes"; then
      if test "x$axes_cv_log4cxx" != "xno"; then
      	 axes_log4espp="log4cxx"
      elif test "x$axes_cv_log4cpp" != "xno"; then
         axes_log4espp="log4cpp"
      else
         axes_log4espp="generic"
      fi
elif test "x$with_log4espp" = "xno"; then
     axes_log4espp="no"
elif test "x$with_log4espp" = "xgeneric"; then
      axes_log4espp="generic"
elif test "x$with_log4espp" = "xlog4cpp"; then
      if test "x$axes_cv_log4cpp" != "xno"; then
      	 axes_log4espp="log4cpp"
      else
         AC_MSG_RESULT([ERROR])
         AC_MSG_ERROR([log4cpp was specified but not found!]);
      fi
elif test "x$with_log4espp" = "xlog4cxx"; then
      if test "x$axes_cv_log4cxx" != "xno"; then
      	 axes_log4espp="log4cxx"
      else
         AC_MSG_RESULT([ERROR])
         AC_MSG_ERROR([log4cxx was specified but not found!]);
      fi
else
      AC_MSG_RESULT([ERROR])
      AC_MSG_ERROR([an unknown logger type was specified as argument to --with-log4espp])
fi

dnl Now evaluate the result

if test "x$axes_log4espp" = "xlog4cxx"; then
     LOG4ESPP_LIBS="-l$axes_cv_log4cxx"
     AC_DEFINE(LOG4ESPP_USE_LOG4CXX,1,[define if log4cxx shall be used by log4espp])
 elif test "x$axes_log4espp" = "xlog4cpp"; then
     LOG4ESPP_LIBS="-l$axes_cv_log4cpp"
     AC_DEFINE(LOG4ESPP_USE_LOG4CPP,1,[define if log4cpp shall be used by log4espp])
elif test "x$axes_log4espp" = "xgeneric"; then
     LOG4ESPP_LIBS=""
     AC_DEFINE(LOG4ESPP_USE_GENERIC,1,[define if the generic logger shall be used by log4espp])
else
     LOG4ESPP_LIBS=""
fi
AC_SUBST([LOG4ESPP_LIBS])

LOG4ESPP_USE=$axes_log4espp
AC_SUBST([LOG4ESPP_USE])
AC_MSG_RESULT([$LOG4ESPP_USE])

])
