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

# we use C++ as main code language
# and libtool for building the shared libraries

cat <<EOF
**********************************************************************
*         Checking for compiler and dynamic library support          *
**********************************************************************
EOF

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

AXES_BOOST_UNIT_TEST_FRAMEWORK
AM_CONDITIONAL([HAVE_UNITTEST], [test "x$axes_cv_boost_lib_unit_test_framework" != "xno"])

cat <<EOF
**********************************************************************
*         Checking for optional external libraries and tools         *
**********************************************************************
EOF

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
		 src/cpp/pmi/mpitransmit/Makefile
                 src/cpp/mpi/Makefile
                 src/cpp/hello/Makefile
                 src/cpp/logging/Makefile
                 test/Makefile 
                 test/external/Makefile 
                 test/external/mpi/Makefile 
                 test/external/logging/Makefile 
                 test/external/boost.mpi/Makefile 
                 test/external/boost.python/Makefile 
                 test/external/boost.python/hello/Makefile 
	         ])

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

if test "x$axes_cv_boost_lib_unit_test_framework" = "xyes"; then
echo " * using the boost unit test framework for unit testing"
fi

if DX_TEST_FEATURE([doc]); then
echo " * Doxygen: doxygen was found"
else
echo " * Doxygen: doxygen was not found"
fi
