#
# SYNOPSIS
#
#   AXES_PYTHON
#
# DESCRIPTION
#
#   This macro checks for the Python headers and libraries.
#
#   It recurses through several python versions (from 2.3 to 2.6 in this
#   version), looking for the library name. If it finds it, it looks to find
#   the header files for this version.
#
#   This macro calls
#      AC_SUBST(PYTHON_CPPFLAGS)
#      AC_SUBST(PYTHON_LIBS)
#
#   And defines:
#
#     HAVE_PYTHON
#
#   And sets the ac cache variables:
#
#     axes_cv_python_lib
#     axes_cv_python_include
#       which can be used to check whether python is available. In this
#       case, both are set to "no".
#
# LAST MODIFICATION
#
#   2008-08-19
#
# COPYLEFT
#
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
AC_ARG_WITH([python],
    AS_HELP_STRING([--with-python@<:@=LIB@:>@],
    [use python (default is yes) - you can specify the library of python
    to use, by default the latest available is used]),
    [
    if test "$withval" = "no"; then
        want_python="no"
    elif test "$withval" = "yes"; then
        want_python="yes"
    else
        want_python="yes"
        axes_python_version="$withval"
    fi
    ],
    [ want_python="yes" ])

if test "x$want_python" != "xno"; then

    dnl save flags before tests
    axes_python_saved_cppflags="$CPPFLAGS"
    axes_python_saved_libs="$LIBS"

    AC_CACHE_CHECK(whether the python library is available,
        axes_cv_python_lib,
    [
        dnl which names to probe for the library
        if test "x$axes_python_version" = "x"; then
            axes_python_possible_libs="python2.5 python2.4 python2.3 python"
        else
            axes_python_possible_libs="$axes_python_version"
        fi
        for axes_python_path_tmp in $axes_python_possible_libs; do
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
                axes_cv_python_lib=$axes_python_path_tmp,
                axes_cv_python_lib=no)
           if test "x$axes_cv_python_lib" != "xno"; then break; fi
        done
    ])

    if test "x$axes_cv_python_lib" != "xno"; then
        AC_CACHE_CHECK(for matching Python.h,
            axes_cv_python_include,
        [
            dnl check possible locations
            for axes_python_path_tmp in /usr/include /usr/local/include \
                    /opt/include /opt/local/include; do
                CPPFLAGS="$axes_python_saved_cppflags -I$axes_python_path_tmp/$axes_cv_python_lib"
                AC_LINK_IFELSE(AC_LANG_PROGRAM([[
                        #include <Python.h>
                	]], [[
                		Py_SetProgramName("test");
    	                return 0;
            	    ]]),
                    axes_cv_python_include="$axes_python_path_tmp/$axes_cv_python_lib",
                    axes_cv_python_include=no)
                if test "x$axes_cv_python_include" != "xno"; then break; fi
            done
        ])
    fi

    dnl restore flags
    CPPFLAGS="$axes_python_saved_cppflags"
    LIBS="$axes_python_saved_libs"

    dnl results
    if test "x$axes_cv_python_include" != "x" && \
       test "x$axes_cv_python_lib" != "x"; then
        PYTHON_CPPFLAGS="-I$axes_cv_python_include"
        AC_SUBST(PYTHON_CPPFLAGS)
        PYTHON_LIBS="-l$axes_cv_python_lib"
        AC_SUBST(PYTHON_LIBS)

        AC_DEFINE(HAVE_PYTHON,1,[define if python is available])
    else
        axes_cv_python_include=no
        axes_cv_python_lib=no
    fi
fi
])
