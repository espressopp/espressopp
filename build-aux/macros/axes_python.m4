#
# SYNOPSIS
#
#   AXES_PYTHON
#
# DESCRIPTION
#
#   This macro checks for the Python headers and libraries.
#
#   If python is not found to work as is, this macro iterates over
#   several python versions (from 2.7 to 2.3), looking for the library.
#   If it is found, the macro tries to find matching header files.
#
#   This macro calls
#
#      AC_SUBST(PYTHON_CPPFLAGS)
#      AC_SUBST(PYTHON_LIBS)
#
#   And defines:
#
#     HAVE_PYTHON
#
#   And sets the ac cache variables:
#
#     axes_cv_python_asis
#     axes_cv_python_lib
#     axes_cv_python_include
#
#   axes_cv_python_asis is set to "yes" if python can be used as is. 
#   In that case, the two other variables are set to "yes".
#   If axes_cv_python_asis is "no", then axes_cv_python_lib and
#   axes_cv_python_include will be set to "no" (if python was not
#   found), or to the required library and preprocessor flags
#   (identical to PYTHON_LIBS and PYTHON_CPPFLAGS).
#
# LAST MODIFICATION
#
#   2009-04-01
#
# COPYLEFT
#
#   Copyright (c) 2009 Olaf Lenz <lenzo@mpip-mainz.mpg.de>
#   Copyright (c) 2008 Axel Arnold <Axel.Arnold@scai.fhg.de>
#
#   based on AX_PYTHON, which is copylefted as
#   
#   Copyright (c) 2008 Michael Tindal
#
#   This program is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation; either version 2 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AXES_PYTHON],
[
# first test, whether compiling a Python program works as is
AC_CACHE_CHECK(
    [whether a python program can be linked as is],
    axes_cv_python_asis,
    AC_LINK_IFELSE(AC_LANG_PROGRAM([[
        #include <Python.h>
        ]], [[
        Py_SetProgramName("test");
        return 0;
        ]]),
      [axes_cv_python_asis=yes],
      [axes_cv_python_asis=no]))

# if not, test the details
AS_IF([test "x$axes_cv_python_asis" = "xno"],
  [ 
    AC_ARG_WITH([python],
        AS_HELP_STRING([--with-python@<:@=LIB@:>@],
          [use python (default is yes) - you can specify the library of python
          to use, by default the latest available is used]),
        AS_IF(
          [test "x$withval" = "xno"],
            [want_python="no"],
          [test "x$withval" = "xyes"]
            [want_python="yes"]
          [
            want_python="yes"
            axes_python_version="$withval"
          ]),
        [ want_python="yes" ])

    AS_IF(
      [test "x$want_python" != "xno"],
        [
          # save flags before tests
          axes_python_saved_cppflags="$CPPFLAGS"
          axes_python_saved_libs="$LIBS"

	  # check for the library
          AC_CACHE_CHECK(
            [whether the python library is available],
            axes_cv_python_lib,
            [
              # which names to probe for the library
              AS_IF(
                [test "x$axes_python_version" = "x"],
                  [axes_python_possible_libs="python2.7 python2.6 python2.5 python2.4 python2.3 python"],
                [axes_python_possible_libs="$axes_python_version"])
              for axes_python_path_tmp in $axes_python_possible_libs
              do
		LIBS="$axes_python_saved_libs -l$axes_python_path_tmp"
                AC_LINK_IFELSE(AC_LANG_PROGRAM([[
                    #ifdef __cplusplus
                    extern "C" {
                    #endif
                      void Py_Main();
                    #ifdef __cplusplus
                    }
                    #endif
                    ]], [[
                    Py_Main();
                    return 0;
                    ]]),
                  [axes_cv_python_lib="$axes_python_path_tmp"],
                  [axes_cv_python_lib=no])
                AS_IF(
                  [test "x$axes_cv_python_lib" != "xno"],
                    [break])
              done
            ])

          # check for matching headers
	  AS_IF([test "x$axes_cv_python_lib" != "xno"],
            [AC_CACHE_CHECK(
              [for matching python headers],
              axes_cv_python_include,
              [
                dnl check include paths from CPPFLAGS
                for axes_python_path_tmp in $CPPFLAGS
                do
                  dnl translate -I flags into paths
                  AS_CASE(["$axes_python_path_tmp"],
                    [-I*], 
                      [axes_python_abs_path=`cd ${axes_python_path_tmp#-I} && pwd`
                       axes_python_try_roots="$axes_python_try_roots $axes_python_abs_path"],
                    []
                  )
                done
               
                axes_python_try_roots="$axes_python_try_roots /usr/include /usr/local/include"

                dnl check possible locations
                for axes_python_path_tmp in $axes_python_try_roots
                do
                  CPPFLAGS="$axes_python_saved_cppflags -I$axes_python_path_tmp/$axes_cv_python_lib"
                  AC_LINK_IFELSE(AC_LANG_PROGRAM([[
                      #include <Python.h>
                      ]], [[
                      Py_SetProgramName("test");
    	              return 0;
            	      ]]),
                    [axes_cv_python_include="$axes_python_path_tmp/$axes_cv_python_lib"],
                    [axes_cv_python_include=no])
                  AS_IF(
                    [test "x$axes_cv_python_include" != "xno"],
                      [break])
                done
              ])])

          # restore flags
          CPPFLAGS="$axes_python_saved_cppflags"
          LIBS="$axes_python_saved_libs"
    ]) # want_python
]) # python_asis

# results
AS_IF([test "x$axes_cv_python_asis" = "xyes"],
        [ AC_DEFINE(HAVE_PYTHON,1,[define if python is available])
          axes_cv_python_include=yes
          axes_cv_python_lib=yes ],
      [test "x$axes_cv_python_include" != "x" &&
          test "x$axes_cv_python_lib" != "x"],
        [ AC_DEFINE(HAVE_PYTHON,1)
          PYTHON_CPPFLAGS="-I$axes_cv_python_include"
          PYTHON_LIBS="-l$axes_cv_python_lib" ],
      [ axes_cv_python_include=no
        axes_cv_python_lib=no])

AC_SUBST(PYTHON_CPPFLAGS)
AC_SUBST(PYTHON_LIBS)

])
