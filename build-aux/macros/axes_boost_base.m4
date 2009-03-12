#
# SYNOPSIS
#
#   AXES_BOOST_BASE([MINIMUM-VERSION])
#
# DESCRIPTION
#
#   Test for the headers of the Boost C++ library of a particular
#   version (or newer).
#
#   If no path to the installed boost library is given, the macro
#   searches the headers in all paths given in LDFLAGS, in some
#   default locations and it evaluates the BOOST_ROOT environment
#   variable for a staged build of boost.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_CPPFLAGS)
#
#   And defines:
#
#     HAVE_BOOST
#
#   And sets the ac cache variable:
#
#     axes_cv_boost
#       which can be used to check whether the boost headers are available
#
#   Note:
#
#   This macro needs to be called _before_ any test for a specific
#   Boost library.  The libraries do not test independently for the
#   existence of the header files.
#
# LAST MODIFICATION
#
#   2009-03-12
#
# COPYLEFT
#
#   Copyright (c) 2009 Olaf Lenz <lenzo@mpip-mainz.mpg.de>
#   Copyright (c) 2008 Axel Arnold <Axel.Arnold@scai.fhg.de>
#
#   Parts of code are
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty provided
#   the copyright notice and this notice are preserved.

AC_DEFUN([AXES_BOOST_BASE],
[
AC_ARG_WITH([boost],
    AS_HELP_STRING([--with-boost@<:@=INCLUDEDIR@:>@],
    [use boost (default is yes) - you can specify the directory where the
    boost headers are located, otherwise configure will search for them]),
    [
    if test "$withval" = "no"; then
        want_boost="no"
    elif test "$withval" = "yes"; then
        want_boost="yes"
    else
        want_boost="yes"
        axes_boost_user_inc_path="$withval"
    fi
    ],
    [ want_boost="yes" ])

AC_ARG_VAR([BOOST_ROOT],
  [directory where boost is installed])

if test "x$want_boost" = "xyes"; then
    dnl determine required boost version

    boost_lib_version_req=ifelse([$1], ,1.20.0,$1)
    boost_lib_version_req_shorten=`expr $boost_lib_version_req : '\([[0-9]]*\.[[0-9]]*\)'`
    boost_lib_version_req_major=`expr $boost_lib_version_req : '\([[0-9]]*\)'`
    boost_lib_version_req_minor=`expr $boost_lib_version_req : '[[0-9]]*\.\([[0-9]]*\)'`
    boost_lib_version_req_sub_minor=`expr $boost_lib_version_req : '[[0-9]]*\.[[0-9]]*\.\([[0-9]]*\)'`
    if test "x$boost_lib_version_req_sub_minor" = "x" ; then
        boost_lib_version_req_sub_minor="0"
    fi
    WANT_BOOST_VERSION=`expr $boost_lib_version_req_major \* 100000 \+ \
        $boost_lib_version_req_minor \* 100 \+ $boost_lib_version_req_sub_minor`

    dnl and look for that version
    AC_CACHE_CHECK([for boost headers with version >= $boost_lib_version_req],
        axes_cv_boost,
        [
        dnl compile the list of paths that we will search for includes,

        if test "x$axes_boost_user_inc_path" != "x"; then
            dnl only try user given include path, if set
            axes_boost_try_incpaths=$axes_boost_user_inc_path
        else
            dnl this will try it without include path, which should be the first option
            axes_boost_try_incpaths="yes"

            dnl then the BOOST_ROOT environment, if set
            if test "x$BOOST_ROOT" != "x"; then
                if test -d "$BOOST_ROOT"; then
                    for axes_boost_path_tmp in "$BOOST_ROOT"/include "$BOOST_ROOT"/include/boost-*; do
                        if test -d "$axes_boost_path_tmp"; then
                            axes_boost_try_incpaths="$axes_boost_try_incpaths $axes_boost_path_tmp"
                        fi
                    done
                fi
            fi

            dnl which paths we try as basis for the include paths of the boost-headers
            axes_boost_try_roots=

            dnl check include paths from CPPFLAGS
            for axes_boost_path_tmp in $CPPFLAGS; do
                dnl translate -I flags into paths
                case "$axes_boost_path_tmp" in
                -I*) axes_boost_abs_path=`cd ${axes_boost_path_tmp#-I} && pwd`
                     axes_boost_try_roots="$axes_boost_try_roots $axes_boost_abs_path" ;;
                esac
            done
            axes_boost_try_roots="$axes_boost_try_roots /usr/include /usr/local/include"

            dnl check include directories and maybe boost-<version> subdirs
            for axes_boost_path_tmp in $axes_boost_try_roots; do
                dnl without version
                if test -d "$axes_boost_path_tmp"; then
                    axes_boost_try_incpaths="$axes_boost_try_incpaths $axes_boost_path_tmp"
                fi
                dnl with version
                for axes_boost_path_tmp_tmp in "$axes_boost_path_tmp/boost-"*; do
                    if test -d "$axes_boost_path_tmp_tmp"; then
                        axes_boost_try_incpaths="$axes_boost_try_incpaths $axes_boost_path_tmp_tmp"
                    fi
                done
            done
        fi
        echo "looking for includes in >$axes_boost_try_incpaths<" >&AS_MESSAGE_LOG_FD

        dnl compile-test the candidates

        dnl save current flags
        axes_boost_cppflags_saved="$CPPFLAGS"

        dnl test include paths for a suitable boost
        succeeded=no
        for axes_boost_path_tmp in $axes_boost_try_incpaths; do
            dnl yes means we get try without an include path
            CPPFLAGS="$axes_boost_cppflags_saved"
            if test $axes_boost_path_tmp != "yes"; then
                CPPFLAGS="-I$axes_boost_path_tmp $CPPFLAGS"
            fi

            AC_LANG_PUSH([C++])
            AC_COMPILE_IFELSE([AC_LANG_PROGRAM([[
               @%:@include <boost/version.hpp>
            ]], [[
                #if BOOST_VERSION >= $WANT_BOOST_VERSION
                // Everything is okay
                #else
                #  error Boost version is too old
                #endif
            ]])],[
                succeeded=yes
            ],[
            ])
            AC_LANG_POP([C++])
            if test $succeeded = yes; then break; fi
        done

        dnl restore flags
        CPPFLAGS="$axes_boost_cppflags_saved"

        dnl and remember the found path
        if test $succeeded = yes; then
            axes_cv_boost=$axes_boost_path_tmp
        else
            axes_cv_boost=no
        fi
    ])

    if test "x$axes_cv_boost" != "xno"; then
        if test "x$axes_cv_boost" != "xyes"; then
            BOOST_CPPFLAGS="-I$axes_cv_boost"
        fi
        AC_SUBST(BOOST_CPPFLAGS)
        AC_DEFINE(HAVE_BOOST,1,[define if the Boost library is available])
    fi
else
    axes_cv_boost="no"
fi
])
