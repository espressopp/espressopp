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

cat >& AS_MESSAGE_FD <<EOF
**********************************************************************
*         Checking for compiler and dynamic library support          *
**********************************************************************
EOF

# we use C++ as main code language
# and libtool for building the shared libraries
AC_PROG_CXX
AC_LANG([C++])

LT_INIT([dlopen])
AM_CONDITIONAL(ENABLE_STATIC, [test "x$enable_static" = "xyes"])
AM_CONDITIONAL(ENABLE_SHARED, [test "x$enable_shared" = "xyes"])

cat >& AS_MESSAGE_FD <<EOF
**********************************************************************
*         Checking system header files                               *
**********************************************************************
EOF

AC_CHECK_HEADERS(sys/time.h sys/resource.h)
AS_IF([test "$ac_cv_header_sys_time_h" != "yes"], 
        [AC_CHECK_HEADER(time.h)])

cat >& AS_MESSAGE_FD <<EOF
**********************************************************************
*       Checking for required external libraries and tools           *
**********************************************************************
EOF

AXES_PYTHON
AS_IF([test "x$axes_cv_python_lib" = "xno" || 
	    test "x$axes_cv_python_include" = "xno"], 
         [AC_MSG_ERROR([could not find python])])

AXES_MPI
AS_IF([test "x$axes_cv_mpi" = "xno"],
       [AC_MSG_ERROR([could not find mpi])])


###################################################################
#                                                                 #
#  mpiexec/mympiexec.sh                                           #
#                                                                 #
###################################################################

AS_IF([test "x$axes_cv_mpi" != "xno"],[
  AC_ARG_VAR([MPIEXEC], [MPI command mpiexec])
  AS_IF([test "x$MPIEXEC" = "x"],
    [AC_PATH_PROG([MPIEXEC], [mpiexec], [no])])

  AC_MSG_CHECKING([for mympiexec])
  AC_ARG_WITH([mympiexec],
     AS_HELP_STRING([--with-mympiexec@<:@=script@:>@],
      [specify the program or script that should be used to run the parallel build tests.
Defaults to BUILDDIR/mympiexec.sh]),
      [AS_IF([test "$withval" = "no"],
         [MYMPIEXEC=""],
	 [MYMPIEXEC="$withval"
          dir=`AS_DIRNAME([$MYMPIEXEC])`
	  AS_IF([test "x$dir" = "x."],
                  [MYMPIEXEC="`pwd`/$MYMPIEXEC"])
         ])],
      [ MYMPIEXEC="`pwd`/mympiexec.sh" ])
  AC_MSG_RESULT($MYMPIEXEC)
  AC_SUBST(MYMPIEXEC)
])


###################################################################
#                                                                 #
#  boost libs
#                                                                 #
###################################################################

AC_MSG_CHECKING([whether to use an external boost library])
AC_ARG_WITH([external-boost],
    AS_HELP_STRING([--with-external-boost=BOOST_MPI_PYTHON_MODULE],
      [specify whether to use external boost libraries (instead of the
       internal ones). If an external boost library is used, specify the
       path of the boost.mpi python module.])
    ,, [with_external_boost="no"])
AS_IF([test "x$with_external_boost" = "xno"],
    AC_MSG_RESULT([no]), AC_MSG_RESULT([yes]))

AC_MSG_CHECKING([where to find the boost.mpi python module])
AS_IF([test "x$with_external_boost" = "xno"],
    [boost_internal_dir=`pwd`"/contrib/boost-am/"
     ESPRESSO_BOOST_LIBS="$boost_internal_dir/libespresso_boost.la"
     BOOST_UNITTEST_LIBS="$boost_internal_dir/libboost_unit_test_framework.la"
     BOOST_MPI_PYTHON_DIR="$boost_internal_dir/.libs"
     BOOST_MPI_PYTHON_LIBS="$boost_internal_dir/libboost_mpi_python.la"
     AC_MSG_RESULT([(internal dir)])
    ],
  [test "x$with_external_boost" = "xyes"],
    [AC_MSG_ERROR([boost.mpi python module was not specified in
     --with-external-boost.])],
  [test ! -e "$with_external_boost"], 
    [AC_MSG_ERROR([file $with_external_boost does not exist! Please specify the complete path of the boost.mpi python module])],
  dnl TODO: Check whether a boost program compiles as is
  [ # else
    BOOST_MPI_PYTHON_DIR=`AS_DIRNAME($with_external_boost)`
    BOOST_MPI_PYTHON_LIBS="$with_external_boost"
    AC_MSG_RESULT([$with_external_boost])
  ]
)

# Test using the unittest framework
AS_IF([test "x$with_external_boost" != "xno"],[
  # test link a program using the dynamic linking
  AC_MSG_CHECKING([whether a boost unittest can be dynamically linked])
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE(AC_LANG_SOURCE([[
        @%:@define BOOST_TEST_DYN_LINK
        @%:@define BOOST_TEST_MODULE conftest
        @%:@include <boost/test/unit_test.hpp>
        BOOST_AUTO_TEST_CASE(conftestcase) {
          BOOST_CHECK(true);
        }
  ]]), 
  axes_boost_unit_test_framework_succeeded=yes,
  axes_boost_unit_test_framework_succeeded=no)
  AC_LANG_POP([C++])
  AC_MSG_RESULT([$axes_boost_unit_test_framework_succeeded])

  # define BOOST_TEST_DYN_LINK if required
  if test "x$axes_boost_unit_test_framework_succeeded" = "xyes"; then
    AC_DEFINE([BOOST_TEST_DYN_LINK], 1, 
    [whether to link the unit test framework library dynamically])
  fi
])

AM_CONDITIONAL([EXTERNAL_BOOST],[test "x$external_boost" = "xyes"])
AC_SUBST(ESPRESSO_BOOST_LIBS)
AC_SUBST(BOOST_UNITTEST_LIBS)
AC_SUBST(BOOST_MPI_PYTHON_LIBS)
AC_SUBST(BOOST_MPI_PYTHON_DIR)

AS_IF([test "x$external_boost" != "xyes"], 
  [AM_CPPFLAGS="$AM_CPPFLAGS -I\$(top_srcdir)/contrib/boost-am/"])

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

AC_ARG_VAR([DOXYGEN], [Doxygen code documentation tool])
AC_PATH_PROG([DOXYGEN],[doxygen],[no])
AM_CONDITIONAL(HAVE_DOXYGEN,[test "x$DOXYGEN" != "xno"])

cat >& AS_MESSAGE_FD <<EOF
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

# determine CPPFLAGS 
AM_CPPFLAGS="-I\$(top_srcdir)/src -I\$(top_srcdir)/src/include \$(PYTHON_CPPFLAGS)"
AS_IF([test "x$external_boost" != "xyes"],
   [
     AM_CPPFLAGS="$AM_CPPFLAGS -I\$(top_srcdir)/contrib/boost-am/"
     UTEXEC_CPPFLAGS="-I\$top_srcdir/contrib/boost-am" 
     UTEXEC_LDFLAGS="-L\$top_builddir/contrib/boost-am"
   ])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([UTEXEC_CPPFLAGS])
AC_SUBST([UTEXEC_LDFLAGS])

################################################################### 
#                                                                 #
#  Generation of output                                           #
#                                                                 #
################################################################### 

AC_CONFIG_HEADERS([src/acconfig.hpp])

AC_CONFIG_FILES([Makefile 
                 contrib/Makefile
                 contrib/boost-am/Makefile
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
                [chmod ugo+x bin/utexec.sh])
AC_CONFIG_FILES([bin/espp_mpiexec.sh],
                [chmod ugo+x bin/espp_mpiexec.sh])
AC_CONFIG_FILES([test/runtests.sh],
                [chmod ugo+x test/runtests.sh])

# espresso setup scripts
AC_CONFIG_FILES([
       espresso_setup.sh:espresso_setup.sh.build.in
       espresso_setup.sh.install:espresso_setup.sh.install.in
       espresso_setup.csh:espresso_setup.csh.build.in
       espresso_setup.py:espresso_setup.py.build.in
       espresso_setup.py.install:espresso_setup.py.install.in
       ])

# create a link to the src dir by the name of espresso, so 
AC_CONFIG_COMMANDS([espresso],
	[AS_IF([test ! -e espresso],[$LN_S src espresso],
               [test ! -h espresso], [
	         AC_MSG_WARN(["espresso exists and is not a link, please remove"])
		 exit 1
               ])])

# now make sure that all config files are generated

AC_OUTPUT

################################################################### 
#                                                                 #
#  Summary                                                        #
#                                                                 #
################################################################### 

test "x$EXTRA_PREFIXES" != "x" || EXTRA_PREFIXES="(none)"
test "x$MPI_LIBS" != "x" || MPI_LIBS="(none required)"
test "x$DOXYGEN" != "xno" || DOXYGEN="(not available)"
AS_IF([test "x$MYMPIEXEC" = "x"],
        [MYMPIEXEC="disabled by user, no parallel tests will be run"],
        [MYMPIEXEC="will try to use $MYMPIEXEC"])
test "x$AM_CPPFLAGS" != "x" || AM_CPPFLAGS="(empty)"

cat >config.summary << EOF


**********************************************************************
*                         Configuration summary                      *
**********************************************************************
*****************
** COMPILATION **
*****************
 * C++-compiler: CXX=$CXX
 * User provided compilation flags:
   - EXTRA_PREFIXES=$EXTRA_PREFIXES
   - CPPFLAGS=$CPPFLAGS
   - AM_CPPFLAGS=$AM_CPPFLAGS (automatically determined by configure)
   - CXXFLAGS=$CXXFLAGS
   - LDFLAGS=$LDFLAGS
 * Library types to build:
   - static libs: $enable_static
   - shared libs: $enable_shared

*****************
** LIBRARIES   **
*****************
 * Python:
   - Headers: $PYTHON_CPPFLAGS
   - Libraries: $PYTHON_LIBS
 * MPI:
   - Libraries: $MPI_LIBS
EOF

AS_IF([test "x$external_boost" = "xyes"], [
cat >>config.summary << EOF
 * Boost: external
   - boost.MPI python dir: $BOOST_MPI_PYTHON_DIR
EOF
],
[
cat >>config.summary << EOF
 * Boost: internal
EOF
])

cat >>config.summary << EOF
 * Log4ESPP: $LOG4ESPP_USE
   - Libraries: $LOG4ESPP_LIBS

*****************
** OTHER TOOLS **
*****************
 * mpiexec:
    - $MYMPIEXEC
    - fallback: $MPIEXEC
 * Doxygen: $DOXYGEN
EOF

cat config.summary >&AS_MESSAGE_FD
cat config.summary >&AS_MESSAGE_LOG_FD

cat >&AS_MESSAGE_FD <<EOF

READY TO BUILD $PACKAGE_STRING

EOF
