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
AC_CONFIG_SRCDIR([src/_espresso.cpp])

# we use automake
AM_INIT_AUTOMAKE([gnu subdir-objects dist-bzip2])
AM_PROG_CC_C_O

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
*         Checking system header files                               *
**********************************************************************
EOF

AC_CHECK_HEADERS(sys/time.h sys/resource.h)
if test "$ac_cv_header_sys_time_h" != "yes"; then
AC_CHECK_HEADER(time.h)
fi

cat <<EOF
**********************************************************************
*       Checking for required external libraries and tools           *
**********************************************************************
EOF
AXES_MPI
if [test "x$axes_cv_mpi" = "xno"]; then
   AC_MSG_ERROR([could not find mpi])   
fi

AXES_PYTHON
if test "x$axes_cv_python_lib" = "xno" || \
   test "x$axes_cv_python_include" = "xno"; then
   AC_MSG_ERROR([could not find python])
fi

# Define only temporary, until all traces of HAVE_MPI and HAVE_PYTHON
# are gone 
AH_BOTTOM([
#define IF_MPI(code) code
#define IF_PYTHON(code) code
])

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

if test "x$axes_cv_mpi" != "xno"; then
   AXES_BOOST_MPI
   if test "x$axes_cv_boost_lib_mpi" = "xno"; then
      AC_MSG_ERROR([could not find the boost mpi library.
Note that this can also mean that the compiler $CXX you are using
to compile ESPResSo++ is not the same you used to build boost.mpi.])
   fi
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


if test "x$axes_cv_mpi" != "xno"; then
  AC_PATH_PROG([MPIEXEC], [mpiexec], [no])
  AC_SUBST(MPIEXEC)

  AC_MSG_CHECKING([for mympiexec])
  AC_ARG_WITH([mympiexec],
     AS_HELP_STRING([--with-mympiexec@<:@=script@:>@],
      [specify the program or script that should be used to run the parallel build tests.
Defaults to BUILDDIR/mympiexec.sh]),
      [
      if test "$withval" = "no"; then
          MYMPIEXEC="" 
      else
          MYMPIEXEC="$withval"

          dir=`AS_DIRNAME([$MYMPIEXEC])`
          if test "x$dir" = "x."; then
            dnl not an absolute path, make absolute
            MYMPIEXEC="`pwd`/$MYMPIEXEC"
          fi
      fi
      ],
      [ MYMPIEXEC="`pwd`/mympiexec.sh" ])
  AC_MSG_RESULT($MYMPIEXEC)
  AC_SUBST(MYMPIEXEC)
fi

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

AC_PATH_PROG(DOXYGEN,doxygen,no)
AM_CONDITIONAL(HAVE_DOXYGEN,[test .$DOXYGEN != .no])

cat <<EOF
**********************************************************************
*                        Writing output files                        *
**********************************************************************
EOF


################################################################### 
#                                                                 #
#  SET ADDITIONAL MAKEFILE VARIABLES                              #
#                                                                 #
################################################################### 
# required by build-aux/macros/axes_python.am
pythondir=${libdir}/${PACKAGE_TARNAME}/espresso
AC_SUBST(pythondir)
# required by src/Makefile.am
AC_SUBST([abs_top_builddir])
AM_CPPFLAGS="-I\$(top_srcdir)/src -I\$(top_srcdir)/src/include $BOOST_CPPFLAGS"
AM_LDFLAGS="$BOOST_LDFLAGS"
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])

################################################################### 
#                                                                 #
#  Generation of output                                           #
#  Makefiles in all subdirectories                                #
#  src/acconfig.hpp                                           #
#                                                                 #
################################################################### 

AC_CONFIG_HEADERS([src/acconfig.hpp])

AC_CONFIG_FILES([Makefile 
                 src/Makefile 
		 src/include/Makefile
                 src/esutil/Makefile
                 src/esutil/unittest/Makefile
                 src/hello/Makefile
                 src/particles/Makefile
                 src/particles/unittest/Makefile
                 src/pairs/Makefile
                 src/pairs/unittest/Makefile
                 src/interaction/Makefile
                 src/bc/Makefile
                 src/bc/unittest/Makefile
                 src/integrator/Makefile
		 doc/Makefile
		 doc/dg/Makefile
                 test/Makefile
                 test/modules/Makefile
                 test/modules/forceloop/Makefile
                 test/external/Makefile 
                 test/external/mpi/Makefile 
                 test/external/logging/Makefile 
                 test/external/boost.mpi/Makefile 
                 test/external/boost.python/Makefile 
                 ])

AC_CONFIG_FILES([bin/utexec.sh],
                [chmod +x bin/utexec.sh])
AC_CONFIG_FILES([bin/espp_mpiexec.sh],
                [chmod +x bin/espp_mpiexec.sh])
AC_CONFIG_FILES([test/runtests.sh],
                [chmod +x test/runtests.sh])

# now make sure that all config files are generated

AC_OUTPUT

cat << EOF


**********************************************************************
*                         Configuration summary                      *
**********************************************************************
*****************
** COMPILATION **
*****************
 * C++-compiler: $CXX
 * Python:
   - Headers: $PYTHON_CPPFLAGS
   - Libraries: $PYTHON_LIBS
 * Boost:
EOF

if test "x$BOOST_CPPFLAGS" != "x"; then
    echo "   - Headers: $BOOST_CPPFLAGS"
else
    echo "   - Headers: default include path"
fi

cat <<EOF
   - Library suffix: ${axes_cv_boost_lib_suffix#-}
   - Linker flags: $BOOST_LDFLAGS
   - boost.Python libraries: $BOOST_PYTHON_LIBS
EOF

dnl Unit test framework
if test "x$BOOST_UNIT_TEST_FRAMEWORK_LIBS" != "x"; then
echo "   - boost.Test libraries: $BOOST_UNIT_TEST_FRAMEWORK_LIBS"
else
echo "   - boost.Test libraries: not available"
fi

cat <<EOF
 * Log4ESPP: $axes_log4espp
   - Libraries: $LOG4ESPP_LIBS
 * Final compilation flags:
   - Preprocessor flags (CPPFLAGS): $CPPFLAGS
   - Compiler flags (CXXFLAGS): $CXXFLAGS
   - Linker flags (LDFLAGS): $LDFLAGS

EOF

cat <<EOF
*****************
** OTHER TOOLS **
*****************
 * mpiexec:
EOF
dnl MPI
if test "$MYMPIEXEC" = ""; then
    echo "   - disabled by user, no parallel tests will be run"
else
    echo "   - will try to use $MYMPIEXEC"
    if ! test -x "$MYMPIEXEC"; then
        echo "   - (which does not exist or is not executable at the moment)"
    fi
    if test "$MPIEXEC" != "no"; then
        echo "   - fallback: $MPIEXEC"
    else
        echo "   - no fallback"
    fi
fi

dnl Doxygen
if test .$DOXYGEN != .no; then
echo " * Doxygen: $DOXYGEN"
else
echo " * Doxygen: doxygen was not found"
fi

echo
echo "READY TO BUILD $PACKAGE_STRING"
echo
