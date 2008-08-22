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

AC_CONFIG_HEADERS(src/cpp/acconfig.h)

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
#  AX_PYTHON sets the following variables                         #
#                                                                 #
#     PYTHON_BIN         : name of the python executable          #
#     PYTHON_INCLUDE_DIR : directory holding the header files     #
#     PYTHON_LIB         : name of the Python library.            #
#                                                                 #
#    Attention: use -I$(PYTHON_INCLUDE_DIR) in Makefile's         #
#                                                                 #
################################################################### 

AX_PYTHON

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

AX_BOOST_BASE(1.35)

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

AC_CONFIG_FILES([Makefile 
                 src/Makefile 
                 src/cpp/Makefile
                 src/cpp/espresso/Makefile 
                 src/cpp/espresso/kernel/Makefile 
                 src/cpp/espresso/kernel/test/Makefile
	         ])

# now make sure that all config files are generated

AC_OUTPUT

