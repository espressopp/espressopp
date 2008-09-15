################################################################### 
#                                                                 #
#  configure.ac for ESPResSo++                                    #
#                                                                 #
#  - use ./autogen.sh to generate the configure script            #
#                                                                 #
################################################################### 

AC_INIT([espressopp], [0.1])
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

AC_PROG_CXX
AC_LANG([C++])

AC_PROG_LIBTOOL

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

AXES_MPI
if test "x$axes_cv_mpi" = "xno"; then
   AC_MSG_ERROR([could not compile an MPI program])
fi

AXES_BOOST_MPI
if test "x$axes_cv_boost_mpi" = "xno"; then
  if test "x$axes_cv_boost_mpi_links" = "xno"; then
      AC_MSG_ERROR([could not link against the boost MPI library.
Are you sure you are using the same MPI compiler/library that was used to
compile the boost mpi library?])
  else
      AC_MSG_ERROR([could not find the boost MPI library])
  fi
fi

AXES_PYTHON
if test "x$axes_cv_python_lib" = "xno" -o \
	"x$axes_cv_python_include" = "xno"; then
   AC_MSG_ERROR([could not find python])
fi

AXES_BOOST_PYTHON
if test "x$axes_cv_boost_lib_python" = "xno"; then
    AC_MSG_ERROR([could not find the boost python library])
fi

AXES_BOOST_UNIT_TEST_FRAMEWORK
if test "x$axes_cv_boost_lib_unit_test_framework" = "xno"; then
    AC_MSG_ERROR([could not find the boost unit test framework])
fi

###################################################################
#                                                                 #
#  check for log4cpp, or log4cxx if log4cpp is not found          #
#                                                                 #
###################################################################

AXES_LOG4CXX
AXES_LOG4CPP

if test "x$axes_cv_log4cxx" != "xno"; then
   axes_log4espp="LOG4CXX"
   AC_SUBST(LOG4ESPP_LIBS,"-l$axes_cv_log4cxx")
   AC_DEFINE(LOG4ESPP_USE_LOG4CXX,1, [define to use log4cxx for logging])
elif test "x$axes_cv_log4cpp" != "xno"; then
   axes_log4espp="LOG4CPP"
   AC_SUBST(LOG4ESPP_LIBS,"-l$axes_cv_log4cpp")
   AC_DEFINE(LOG4ESPP_USE_LOG4CPP,1, [define to use log4cpp for logging])
else
   axes_log4espp="generic logger"
   AC_SUBST(LOG4ESPP_LIBS,"")
fi

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
                 src/cpp/espresso/Makefile 
                 src/cpp/espresso/kernel/Makefile 
                 src/cpp/espresso/kernel/test/Makefile
                 examples/Makefile 
                 examples/cpp/Makefile
	         ])

# now make sure that all config files are generated

AC_OUTPUT

cat << EOF


************************
* INSTALLATION SUMMARY

EOF

echo "* BOOST:"
if test "x$axes_cv_boost" != "xyes"; then
    echo "  - headers located at $axes_cv_boost"
else
    echo "  - headers located in the default includes"
fi
echo "  - flavor ${axes_cv_boost_lib_suffix#-}"
echo

echo "* LOG4ESSP (logging for Espresso): uses $axes_log4espp"
echo

if DX_TEST_FEATURE([doc]); then
echo "* Documentation: uses doxygen"
echo
else
echo "* Documentation: not built"
echo
fi

cat << EOF
*
************************
EOF
