################################################################### 
#                                                                 #
#  configure.ac                                                   #
#                                                                 #
#  Input file for autoconf + automake                             #
#                                                                 #
#  - use script autogen.sh for the initial configuration          #
#                                                                 #
################################################################### 

AC_INIT(src/cpp/espresso/kernel/BasicProperty.cpp)
AM_INIT_AUTOMAKE(EspressoPP,0.1)
AM_CONFIG_HEADER(config.h)

################################################################### 
#                                                                 #
#  MPI compiler driver for C++                                    #
#                                                                 #
################################################################### 

AC_LANG([C++])
ACX_MPI
CXX="$MPICXX"

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

AM_PATH_PYTHON([2.2])

################################################################### 
#                                                                 #
#  BOOST + PYTHON                                                 #
#                                                                 #
################################################################### 

AX_PYTHON
AX_BOOST_BASE(1.35)
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS -I$PYTHON_INCLUDE_DIR"


################################################################### 
#                                                                 #
#  BOOST UNIT TEST FRAMEWORK (AX_BOOST_UNIT_TEST_FRAMEWORK)       #
#                                                                 #
#    - test for Unit_Test_Framework library                       #
#    - sets HAVE_BOOST_UNIT_TEST_FRAMEWORK                        #
#    - calls AC_SUBST(BOOST_UNIT_TEST_FRAMEWORK_LIB)              #
#                                                                 #
################################################################### 

AX_BOOST_UNIT_TEST_FRAMEWORK

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

AC_OUTPUT(Makefile 
          src/Makefile 
          src/cpp/Makefile
          src/cpp/espresso/Makefile 
          src/cpp/espresso/kernel/Makefile 
          src/cpp/espresso/kernel/test/Makefile
	  )
