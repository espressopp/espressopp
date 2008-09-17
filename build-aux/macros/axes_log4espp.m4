# DESCRIPTION
#   
#   Chooses the logger library for LOG4ESPP.
#
#   This macro calls:
#
#     AC_SUBST(LOG4ESPP_LIBS)
#
#   And sets the ac variables:
#
#     axes_log4cxx
#
# LAST MODIFICATION
#
#   2008-09-17
#
# COPYLEFT
#
#   Copyright (c) 2008 Olaf Lenz <lenzo@mpip-mainz.mpg.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_LOG4ESPP],
[
AC_REQUIRE([AXES_LOG4CXX])
AC_REQUIRE([AXES_LOG4CPP])

AC_ARG_WITH([log4espp],
   AS_HELP_STRING([--with-log4espp=LOGGER], 
[use the specified logger library (default: auto). 
Possible values are "log4cxx", "log4cpp", "generic" and "auto".])
  ,[log4espp_wanted=$withval], [log4espp_wanted="auto"])

AC_MSG_CHECKING([which logger library to use])

dnl Test which logger to take
if test "x$log4espp_wanted" = "xauto" || test "x$log4espp_wanted" = "xyes"; then
      if test "x$axes_cv_log4cxx" != "xno"; then
      	 axes_log4espp="log4cxx"
      elif test "x$axes_cv_log4cpp" != "xno"; then
         axes_log4espp="log4cpp"
      else
         axes_log4espp="generic"
      fi
elif test "x$log4espp_wanted" = "xgeneric" || test "x$log4espp_wanted" = "xno"; then
      axes_log4espp="generic"
elif test "x$log4espp_wanted" = "xlog4cpp"; then
      if test "x$axes_cv_log4cpp" != "xno"; then
      	 axes_log4espp="log4cpp"
      else
         AC_MSG_RESULT([ERROR])
         AC_MSG_ERROR([log4cpp was specified but not found!]);
      fi
elif test "x$log4espp_wanted" = "xlog4cxx"; then
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
else
     LOG4ESPP_LIBS=""
fi

AC_MSG_RESULT($axes_log4espp)

AC_SUBST(LOG4ESPP_LIBS)
])
