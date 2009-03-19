#
# SYNOPSIS
#
#   AXES_BOOST_LIB(libname, LIBNAME, program)
#
# DESCRIPTION
#
#   Test for a library from the Boost C++ libraries. The parameter program
#   is the program code used to test compilation / linking. depend_libs
#   specifies 
#   The macro requires a preceding call to AXES_BOOST_BASE. Typically,
#   you do not want to use AXES_BOOST_LIB directly, but rather one of
#   the AXES_BOOST_SERIALIZATION, AXES_BOOST_SIGNAL etc macros, which
#   provide the test program code and additional tests, if necessary.
#
#   This macro calls:
#
#     AC_SUBST(BOOST_<libname>_LIBS)
#     AC_SUBST(BOOST_LDFLAGS)
#
#   And defines:
#
#     HAVE_BOOST_<libname>
#
#   And sets the ac variables:
#
#     axes_cv_boost_<libname>
#     axes_cv_boost_lib_suffix
#     axes_cv_boost_lib_path
#
#   Note: <libname> and <LIBNAME> must be literals, not variables. They should
#     be identical, except for capitalization.
#
# LAST MODIFICATION
#
#   2009-03-11
#
# COPYLEFT
#
#   Copyright (c) 2009 Olaf Lenz <lenzo@mpip-mainz.mpg.de>
#   Copyright (c) 2008 Axel Arnold <axel.arnold@scai.fraunhofer.de>
#
#   based on AX_BOOST_SERIALIZATION, copylefted as follows:
#
#   Copyright (c) 2008 Thomas Porschberg <thomas@randspringer.de>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AXES_BOOST_LIB],
[
AC_REQUIRE([AXES_BOOST_BASE])
AC_ARG_WITH([boost-$1],
AS_HELP_STRING([--with-boost-$1@<:@=special-lib@:>@],
               [use the $1 library from boost - it is possible to specify a certain library for the linker
                    e.g. --with-boost-$1=boost_$1-gcc-mt-d-1_35_0 ]),
    [
    if test "$withval" = "no"; then
        want_boost="no"
        axes_boost_user_lib=""
    elif test "$withval" = "yes"; then
        want_boost="yes"
        axes_boost_user_lib=""
    else
        want_boost="yes"
        axes_boost_user_lib="$withval"
    fi
    ],
    [want_boost="yes"; axes_boost_user_lib=""]
)

AC_ARG_WITH([boost-flavor],
    AS_HELP_STRING([--with-boost-flavor=FLAVOR],
    [specify the flavor of boost libraries to use, that is the name extension,
     e.g. 'gcc-mt-1_35'; otherwise, the flavor is guessed when testing for the first library]),
    [
    if test "$withval" = "no"; then
        axes_boost_user_lib_suffix="none"
    elif test "$withval" = ""; then
        axes_boost_user_lib_suffix="none"
    else
        axes_boost_user_lib_suffix="-$withval"
    fi
    ], [])

dnl make proper variables from the cache ones
AS_VAR_PUSHDEF([cv_boost_lib], [axes_cv_boost_lib_$1])

if test "x$axes_cv_boost" != "xno" && test "x$want_boost" != "xno" ; then
    AC_CACHE_CHECK(whether the boost $1 library is available,
                   cv_boost_lib,
        [
        dnl step 0: try to guess possible paths, if not specified

        dnl compile the list of paths that we will search for the library
        dnl Also, we base the name on the guessed include directory -
        dnl again, mixing headers and libs of different versions
        dnl probably won't work

        if test "x$axes_cv_boost_lib_path" != "x"; then
            dnl the path was already found, we just use that again
            dnl - auto guessing should really not mix boost-installations.
            axes_boost_try_libpaths=$axes_cv_boost_lib_path
        else
            dnl first try without a path
            axes_boost_try_libpaths="yes"

            if test "x$axes_cv_boost" = "xyes"; then
                dnl the headers are in one of the user-specified includes,
                dnl so check, whether we can get from the user-specified
                dnl lib paths or the compiler default paths

                dnl check lib paths from LDFLAGS
                for axes_boost_path_tmp in $LDFLAGS; do
                    dnl translate -L flags into paths
                    case "$p" in
                    -L*) axes_boost_abs_path=`cd ${axes_boost_path_tmp#-L} && pwd`
                         axes_boost_try_libpaths="$axes_boost_try_libpaths $axes_boost_abs_path" ;;
                    esac
                done

                dnl append compiler paths
                axes_boost_try_libpaths="$axes_boost_try_libpaths /usr/lib /usr/local/lib"
            else
                dnl otherwise, we have some path, so only use paths derived from
                dnl that
                axes_boost_try_roots="\
                    `echo $axes_cv_boost | sed 's,/include$,,'` \
                    `echo $axes_cv_boost | sed 's,/include/boost-[[^/]]*$,,'` \
                    $axes_cv_boost/stage"
                dnl but here, we need to guess the library subdirectory
                for axes_boost_root_tmp in $axes_boost_try_roots; do
                    for axes_boost_pathsuffix in lib lib64 lib32 shlib; do
                        axes_boost_path_tmp="$axes_boost_root_tmp/$axes_boost_pathsuffix"
                        if test -d "$axes_boost_path_tmp"; then
                            axes_boost_try_libpaths="$axes_boost_try_libpaths\
                                $axes_boost_path_tmp"
                        fi
                    done
               done
            fi
        fi
        echo "looking for lib in >$axes_boost_try_libpaths<" >&AS_MESSAGE_LOG_FD

        dnl step 1: try to guess possible suffixes, if not specified

        if test "x$axes_boost_user_lib" != "x"; then
            dnl enforce use of user-specified library name
            axes_boost_try_suffixes=".user."
    	elif test "x$axes_boost_user_lib_suffix" != "x"; then
            dnl user defined suffix -> use only that, no guessing
            if test $axes_boost_user_lib_suffix != "none"; then
                axes_boost_user_lib="boost_$1$axes_boost_user_lib_suffix"
            else
                axes_boost_user_lib="boost_$1"
            fi
            axes_boost_try_suffixes=".user."
    	elif test "x$axes_cv_boost_lib_suffix" != "x"; then
            dnl already automatically determined suffix -> use only that, no guessing
            if test $axes_cv_boost_lib_suffix != "none"; then
                axes_boost_user_lib="boost_$1$axes_cv_boost_lib_suffix"
            else
                axes_boost_user_lib="boost_$1"
            fi
            axes_boost_try_suffixes=".user."
        else
            axes_boost_try_suffixes=""
            for axes_boost_path_tmp in $axes_boost_try_libpaths; do
	    	echo "Looking for lib in $axes_boost_path_tmp..."
                for axes_boost_lib_tmp in $axes_boost_path_tmp/*boost_$1*; do
		    echo "Found $axes_boost_lib_tmp..."
                    if test ! -f $axes_boost_lib_tmp; then continue; fi
                    axes_boost_suffix_tmp=`echo $axes_boost_lib_tmp | \
                        sed 's,^.*boost_$1,,' | sed 's,[[.]].*$,,'`
                    dnl it is possible that boost is not versioned at all, protect the empty suffix
                    if test "x$axes_boost_suffix_tmp" = "x"; then
                        axes_boost_suffix_tmp="none"
                    fi
                    dnl check that this suffix was not yet specified
                    case " $axes_boost_try_suffixes " in
                        *" $axes_boost_suffix_tmp "*) ;;
                        *) axes_boost_try_suffixes="$axes_boost_try_suffixes $axes_boost_suffix_tmp" ;;
                    esac
                done
            done
        fi
        echo "trying suffixes >$axes_boost_try_suffixes<" >&AS_MESSAGE_LOG_FD

        dnl step 2: compile+link-test the candidate directories and suffixes

        dnl save current flags
        axes_boost_cppflags_saved="$CPPFLAGS"
        axes_boost_ldflags_saved="$LDFLAGS"
        axes_boost_libs_saved="$LIBS"

        dnl the CPPFLAGS for boost we know already
        CPPFLAGS="$BOOST_CPPFLAGS $CPPFLAGS"

        succeeded=no
        for axes_boost_path_tmp in $axes_boost_try_libpaths; do
            dnl yes actually just means we get along without an include path
            LDFLAGS="$axes_boost_ldflags_saved"
            if test $axes_boost_path_tmp != "yes"; then
                LDFLAGS="-L$axes_boost_path_tmp $LDFLAGS"
            fi
            for axes_boost_suffix_tmp in $axes_boost_try_suffixes; do
                dnl .user. means to stick to user-defined name or
                dnl the suffix from previous library tests
                if test $axes_boost_suffix_tmp = ".user."; then
                    axes_boost_lib=$axes_boost_user_lib
                else
                    dnl none tries without suffix
                    if test $axes_boost_suffix_tmp != "none"; then
                        axes_boost_lib="boost_$1$axes_boost_suffix_tmp"
                    else
                        axes_boost_lib="boost_$1"
                    fi
                fi
                LIBS="$axes_boost_libs_saved -l$axes_boost_lib"

                AC_LANG_PUSH([C++])
                AC_LINK_IFELSE($3, succeeded=yes, succeeded=no)
                AC_LANG_POP([C++])

                if test "x$succeeded" = "xyes"; then break; fi
            done
            if test "x$succeeded" = "xyes"; then break; fi
        done

        if test "x$succeeded" = "xyes"; then
            dnl cache the lib name, also for pretty print
            cv_boost_lib="$axes_boost_lib"
        else
            cv_boost_lib=no
        fi

        dnl restore flags
        CPPFLAGS="$axes_boost_cppflags_saved"
        LDFLAGS="$axes_boost_ldflags_saved"
        LIBS="$axes_boost_libs_saved"
    ])


    if test "x$cv_boost_lib" != "xno"; then
        dnl load path for future checks, if not yet set
        dnl if we use a (cached) library, we
        dnl want to use the corresponding suffix and path
        dnl also for all subsequent libraries, again to avoid
        dnl mixing of different boost-versions without explicit
        dnl user request
        dnl if the cache values are not cached, cv_boost_lib is not
        dnl cached, and therefore the _tmp-paths set
        if test "x$axes_cv_boost_lib_path" = "x"; then
            AC_MSG_CHECKING([location of boost libraries])
            AC_CACHE_VAL(axes_cv_boost_lib_path,
                [ axes_cv_boost_lib_path=$axes_boost_path_tmp ])
            if test "x$axes_cv_boost_lib_path" = "xyes"; then
                AC_MSG_RESULT(standard library paths)
            else
                AC_MSG_RESULT($axes_cv_boost_lib_path)
            fi
        fi
        dnl same for suffix
        if test "x$axes_cv_boost_lib_suffix" = "x"; then
            AC_MSG_CHECKING([suffix of boost libraries])
            AC_CACHE_VAL(axes_cv_boost_lib_suffix, [
                dnl if user specified the lib or suffix, use that
                if test "x$axes_boost_suffix_tmp" = "x.user."; then
                    axes_cv_boost_lib_suffix=`echo $axes_boost_user_lib | sed -e 's,^boost_$1,,'`
                else
                    axes_cv_boost_lib_suffix=$axes_boost_suffix_tmp
                fi
            ])
            if test "x$axes_cv_boost_lib_suffix" = "xnone"; then
                AC_MSG_RESULT([unversioned layout])
                AC_MSG_WARN([
******************************************************************************
    Your boost libraries are not versioned. This means that we cannot check if
    the library matches the header files. If you encounter problems with the
    boost libraries, please check manually that your header file version
    matches your library version.
******************************************************************************
                ])
            else
                AC_MSG_RESULT($axes_cv_boost_lib_suffix)
            fi
	fi

        dnl set macro with library name
        BOOST_$2_LIBS="-l$cv_boost_lib"
        AC_SUBST(BOOST_$2_LIBS)
        
    	AC_DEFINE(HAVE_BOOST_$2,1,[define if boost $1 library is available])

        dnl set the LDFLAGS to load the lib path if a non-standard path is needed
        if test "x$BOOST_LDFLAGS" = "x" && test "x$axes_cv_boost_lib_path" != "xyes"; then        
            BOOST_LDFLAGS="-L$axes_cv_boost_lib_path"
            AC_SUBST(BOOST_LDFLAGS)
        fi
    fi
else
    dnl we even don't have boost headers, so also no library
    cv_boost_lib="no"
    AC_MSG_WARN([skipping boost_$1, since boost headers are missing])
fi

AS_VAR_POPDEF([cv_boost_lib])

])
