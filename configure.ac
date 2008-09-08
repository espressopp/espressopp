################################################################### 
#                                                                 #
#  configure.ac                                                   #
#                                                                 #
#  Input file for autoconf + automake                             #
#                                                                 #
#  - use script autogen.sh for the initial configuration          #
#                                                                 #
################################################################### 

AC_INIT(espressopp, 0.1)
AM_INIT_AUTOMAKE()


################################################################### 
#                                                                 #
#  - autoconfig generates the following file with                 #
#    a lot of variables of HAVE_feature                           #
#                                                                 #
#  - this might be included where necessary                       #
#                                                                 #
################################################################### 

AC_CONFIG_HEADERS(src/hpp/acconfig.hpp)

################################################################### 
#                                                                 #
#  MPI compiler driver for C++                                    #
#                                                                 #
################################################################### 

AC_PROG_CXX
AC_LANG([C++])

AC_PROG_INSTALL
AC_PROG_LIBTOOL

################################################################### 
#                                                                 #
#  PYTHON                                                         #
#                                                                 #
#    - automake provides support for Python compilation           #
#    - ACTION-IF-NOT-FOUND defaults to abort                      #
#                                                                 #
################################################################### 

# AM_PATH_PYTHON([2.2])

################################################################### 
#                                                                 #
#  AX_PYTHON sets the following variables                         #
#                                                                 #
#     PYTHON_BIN         : name of the python executable          #
#     PYTHON_INCLUDE_DIR : directory holding the header files     #
#     PYTHON_LIB         : name of the Python library.            #
#                                                                 #
#    Attention: use -I$(PYTHON_INCLUDE_DIR) in Makefile's         #
#                                                                 #
################################################################### 

################################################################### 
#                                                                 #
#  AX_BOOST_BASE sets the following variables                     #
#                                                                 #
#     BOOST_CPPFLAGS :  flags used for compilation                #
#     BOOST_LDFLAGS  :  flags used for linking boost              #
#                                                                 #
#     Note: BOOST_LDFLAGS contains the link path but not          #
#           the libraries itself (that will follow)               #
#                                                                 #
################################################################### 

AXES_BOOST_BASE(1.35)

################################################################### 
#                                                                 #
#  BOOST UNIT TEST FRAMEWORK (AX_BOOST_UNIT_TEST_FRAMEWORK)       #
#                                                                 #
#    - test for Unit_Test_Framework library                       #
#    - sets HAVE_BOOST_UNIT_TEST_FRAMEWORK                        #
#    - sets BOOST_UNIT_TEST_FRAMEWORK_LIB to the needed lib       #
#                                                                 #
################################################################### 

AXES_BOOST_PYTHON
AXES_BOOST_UNIT_TEST_FRAMEWORK

###################################################################
#                                                                 #
#  LOG4CPP                                                        #
#                                                                 #
#    - if found we take it for LOG4ESPP                           #
#                                                                 #
###################################################################

AX_LOG4CPP()

if test "$HAVE_LOG4CPP" == "1"; then
   echo "LOG4CPP is used for LOG4ESSP (logging for Espresso)"
   AC_SUBST(LOG4ESPP_CPPFLAGS,${LOG4CPP_CPPFLAGS})
   AC_SUBST(LOG4ESPP_LDFLAGS,${LOG4CPP_LDFLAGS})
   AC_SUBST(LOG4ESPP_LIBS,"-llog4cpp")
fi

################################################################### 
#                                                                 #
#  LOG4CXX                                                        #
#                                                                 #
#    - if found we take it for LOG4ESPP                           #
#                                                                 #
################################################################### 

AX_LOG4CXX()

if test "$HAVE_LOG4CXX" == "1"; then
   echo "LOG4CXX is used for LOG4ESSP (logging for Espresso)"
   AC_SUBST(LOG4ESPP_CPPFLAGS,${LOG4CXX_CPPFLAGS})
   AC_SUBST(LOG4ESPP_LDFLAGS,${LOG4CXX_LDFLAGS})
   AC_SUBST(LOG4ESPP_LIBS,"-llog4cxx")
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
#  Generation of Makefile in all subdirectories                   #
#                                                                 #
################################################################### 

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

