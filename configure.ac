################################################################### 
#                                                                 #
#  configure.ac for ESPResSo++                                    #
#                                                                 #
#  - use ./autogen.sh to generate the configure script            #
#                                                                 #
################################################################### 

AC_INIT([ESPResSo++], [0.1], [espressopp-user@www.espresso-pp.de], [espressopp])
# basic autoconf configuration
AC_CONFIG_MACRO_DIR(build-aux/macros)
AC_CONFIG_AUX_DIR(build-aux)
AC_CONFIG_SRCDIR([src/cpp/hello.hpp])

# we use automake
AM_INIT_AUTOMAKE

# import paths from the EXTRA_PREFIXES precious variable
AXES_IMPORT_EXTRA_PREFIXES

cat <<EOF
**********************************************************************
*         Checking for compiler and dynamic library support          *
**********************************************************************
EOF

# we use C++ as main code language
# and libtool for building the shared libraries
AC_PROG_CXX
AC_LANG([C++])

AC_PROG_LIBTOOL


cat <<EOF
**********************************************************************
*       Checking for required external libraries and tools           *
**********************************************************************
EOF
AXES_MPI
if test "x$axes_cv_mpi" = "xno"; then
    AC_MSG_ERROR([could not compile an MPI program])
fi

AXES_PYTHON
if test "x$axes_cv_python_lib" = "xno" || \
   test "x$axes_cv_python_include" = "xno"; then
   AC_MSG_ERROR([could not find python])
fi

################################################################### 
#                                                                 #
#  check for various boost libraries. All are required, and       #
#  we bail out if one or more are not found.                      #
#                                                                 #
################################################################### 

AXES_BOOST_BASE([1.35])

if test "x$axes_cv_boost" = "xno"; then
    AC_MSG_ERROR([could not find the boost headers])
fi

AXES_BOOST_SERIALIZATION
if test "x$axes_cv_boost_lib_serialization" = "xno"; then
    AC_MSG_ERROR([could not find the boost serialization library])
fi

AXES_BOOST_MPI
if test "x$axes_cv_boost_lib_mpi" = "xno"; then
    AC_MSG_ERROR([could not find the boost mpi library.
Note that this can also mean that the compiler $CXX you are using
to compile ESPResSo++ is not the same you used to build boost.mpi.])
fi

AXES_BOOST_PYTHON
if test "x$axes_cv_boost_lib_python" = "xno"; then
    AC_MSG_ERROR([could not find the boost python library])
fi

cat <<EOF
**********************************************************************
*         Checking for optional external libraries and tools         *
**********************************************************************
EOF

###################################################################
#                                                                 #
#  mpiexec/mympiexec.sh                                           #
#                                                                 #
###################################################################

AC_PATH_PROG([MPIEXEC], [mpiexec], [no])
AC_SUBST(MPIEXEC)

AC_ARG_WITH([mympiexec],
    AS_HELP_STRING([--with-mympiexec@<:@=script@:>@],
    [specify the program or script that should be used to run the parallel build tests.
The program should understand the basic syntax of the MPI-2 mpiexec command.
Defaults to <builddir>/mympiexec.sh]),
    [
    if test "$withval" = "no"; then
        MYMPIEXEC=""
    else
        MYMPIEXEC="$withval"

        if test "${MYMPIEXEC#/}" = "$MYMPIEXEC"; then
            dnl not an absolute path, make absolute
            MYMPIEXEC="`pwd`/$MYMPIEXEC"
        fi
    fi
    ],
    [ MYMPIEXEC="$ac_abs_top_builddir/mympiexec.sh" ])
AC_SUBST(MYMPIEXEC)

###################################################################
#                                                                 #
#  optional boost packages                                        #
#                                                                 #
###################################################################

AXES_BOOST_UNIT_TEST_FRAMEWORK
AM_CONDITIONAL([HAVE_UNITTEST], [test "x$axes_cv_boost_lib_unit_test_framework" != "xno"])

###################################################################
#                                                                 #
#  check for log4cpp, or log4cxx if log4cpp is not found          #
#                                                                 #
###################################################################

AXES_LOG4ESPP

################################################################### 
#                                                                 #
#  DOXYGEN                                                        #
#                                                                 #
################################################################### 

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(Espresso, doxygen.cfg, doc)

cat <<EOF
**********************************************************************
*                        Writing output files                        *
**********************************************************************
EOF

CPPFLAGS="$CPPFLAGS -I\$(top_srcdir)/src/cpp"

################################################################### 
#                                                                 #
#  Generation of output                                           #
#  Makefiles in all subdirectories                                #
#  src/cpp/acconfig.hpp                                           #
#                                                                 #
################################################################### 

AC_CONFIG_HEADERS([src/cpp/acconfig.hpp])

AC_CONFIG_FILES([Makefile 
                 src/Makefile 
                 src/cpp/Makefile
                 src/cpp/pmi/Makefile
                 src/cpp/mpi/Makefile
                 src/cpp/hello/Makefile
                 src/cpp/logging/Makefile
                 src/python/Makefile
                 test/Makefile
                 test/modules/Makefile
                 test/external/Makefile 
                 test/external/mpi/Makefile 
                 test/external/logging/Makefile 
                 test/external/boost.mpi/Makefile 
                 test/external/boost.python/Makefile 
	         ])

AC_CONFIG_FILES([ bin/espp_mpiexec.sh ],
                [ chmod +x bin/espp_mpiexec.sh ])

# now make sure that all config files are generated

AC_OUTPUT

cat << EOF


**********************************************************************
*                         Configuration summary                      *
**********************************************************************
 * C++-compiler: $CXX
 * Python:
    - headers found in $axes_cv_python_include
    - headers found in $axes_cv_python_lib
 * Boost:
EOF
if test "x$axes_cv_boost" != "xyes"; then
    echo "   - headers found in $axes_cv_boost"
else
    echo "   - headers found in the default includes"
fi
echo "   - flavor ${axes_cv_boost_lib_suffix#-}"

echo " * LOG4ESSP (logging for ESPResSo++): uses $axes_log4espp"

if test "$MYMPIEXEC" = ""; then
echo " * mpiexec disabled by user, no parallel tests will be run"
else
    if test -x "$MYMPIEXEC"; then
echo " * mpiexec: using mympiexec = $MYMPIEXEC"
    else
echo " * mpiexec: using mympiexec = $MYMPIEXEC (not present or not executable at the moment)"
    fi
    if test "$MPIEXEC" != "no"; then
echo "            fallback is $MPIEXEC, if mympiexec is not present"
    else
echo "            no fallback, will not run parallel test if mympiexec is not present"
    fi
fi

if test "x$axes_cv_boost_lib_unit_test_framework" != "xno"; then
echo " * using the boost unit test framework for unit testing"
fi

if DX_TEST_FEATURE([doc]); then
echo " * Doxygen: doxygen was found"
else
echo " * Doxygen: doxygen was not found"
fi
