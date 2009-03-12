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

cat >& AS_MESSAGE_FD <<EOF
**********************************************************************
*         Checking system header files                               *
**********************************************************************])
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
#  check for various boost libraries. All are required, and       #
#  we bail out if one or more are not found.                      #
#                                                                 #
################################################################### 

AXES_BOOST_BASE([1.35])

AS_IF([test "x$axes_cv_boost" = "xno"],
    [AC_MSG_ERROR([could not find the boost headers])])

AXES_BOOST_SERIALIZATION
AS_IF([test "x$axes_cv_boost_lib_serialization" = "xno"],
    [AC_MSG_ERROR([could not find the boost serialization library])])

AXES_BOOST_PYTHON
AS_IF([test "x$axes_cv_boost_lib_python" = "xno"],
    [AC_MSG_ERROR([could not find the boost python library])])

AXES_BOOST_MPI
AS_IF([test "x$axes_cv_boost_lib_mpi" = "xno"],
    [AC_MSG_ERROR([could not find the boost mpi library.
Note that this can also mean that the compiler $CXX you are using
to compile ESPResSo++ is not the same you used to build boost.mpi.])])

cat >& AS_MESSAGE_FD <<EOF
**********************************************************************
*         Checking for optional external libraries and tools         *
**********************************************************************
EOF

###################################################################
#                                                                 #
#  mpiexec/mympiexec.sh                                           #
#                                                                 #
###################################################################

AS_IF([test "x$axes_cv_mpi" != "xno"],[
  AC_MSG_CHECKING([for mpiexec])
  AC_ARG_VAR([MPIEXEC], [MPI command mpiexec])
  AS_IF([test "x$MPIEXEC" = "x"],
    [AC_PATH_PROG([MPIEXEC], [mpiexec], [no])])
  AC_MSG_RESULT($MPIEXEC)

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
AM_CPPFLAGS="-I\$(top_srcdir)/src -I\$(top_srcdir)/src/include $BOOST_CPPFLAGS"
AM_LDFLAGS="$BOOST_LDFLAGS"
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])

################################################################### 
#                                                                 #
#  Generation of output                                           #
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
                [chmod ugo+x bin/utexec.sh])
AC_CONFIG_FILES([bin/espp_mpiexec.sh],
                [chmod ugo+x bin/espp_mpiexec.sh])
AC_CONFIG_FILES([test/runtests.sh],
                [chmod ugo+x test/runtests.sh])

# espresso setup scripts
AC_CONFIG_FILES([
       espresso_setup.sh:espresso_setup.sh.build.in
       espresso_setup.sh.install:espresso_setup.sh.install.in
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
test "x$BOOST_ROOT" != "x" || BOOST_ROOT="(empty)"
test "x$BOOST_CPPFLAGS" != "x" || BOOST_CPPFLAGS="(default include path)"
test "x$BOOST_LDFLAGS" != "x" || BOOST_LDFLAGS="(none)"
test "x$BOOST_UNIT_TEST_FRAMEWORK_LIBS" != "x" || BOOST_UNIT_TEST_FRAMEWORK_LIBS="(not available)"
test "x$BOOST_MPI_PYTHON_DIR" != "x" || BOOST_MPI_PYTHON_DIR="(not set)"
test "x$DOXYGEN" != "xno" || DOXYGEN="(not available)"
AS_IF([test "x$MYMPIEXEC" = "x"],
        [MYMPIEXEC="disabled by user, no parallel tests will be run"],
        [MYMPIEXEC="will try to use $MYMPIEXEC"])
test "x$AM_CPPFLAGS" != "x" || AM_CPPFLAGS="(empty)"
test "x$AM_LDFLAGS" != "x" || AM_LDFLAGS="(empty)"

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
   - CXXFLAGS=$CXXFLAGS
   - LDFLAGS=$LDFLAGS
 * Compilation flags determined by configure:
   - AM_CPPFLAGS=$AM_CPPFLAGS
   - AM_LDFLAGS=$AM_LDFLAGS

*****************
** LIBRARIES   **
*****************
 * Python:
   - Headers: $PYTHON_CPPFLAGS
   - Libraries: $PYTHON_LIBS
 * MPI:
   - Libraries: MPI_LIBS=$MPI_LIBS
 * Boost:
   - BOOST_ROOT=$BOOST_ROOT
   - Headers: $BOOST_CPPFLAGS
   - Library suffix: ${axes_cv_boost_lib_suffix#-}
   - Linker flags: $BOOST_LDFLAGS
   - boost.Python libraries: $BOOST_PYTHON_LIBS
   - boost.MPI libraries: $BOOST_MPI_LIBS
   - boost.Test libraries: $BOOST_UNIT_TEST_FRAMEWORK_LIBS
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
