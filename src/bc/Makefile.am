subpkgincludedir=$(pkgincludedir)/bc

nobase_subpkginclude_HEADERS = \
	BC.hpp \
	PBC.hpp

noinst_LTLIBRARIES = libbc.la

libbc_la_SOURCES = \
	PBC.cpp

check_PROGRAMS = pbcTest

TESTS = $(check_PROGRAMS)

pbcTest_SOURCES = unittest/pbcTest.cpp PBC.cpp

pbcTest_CPPFLAGS = $(BOOST_CPPFLAGS)
pbcTest_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

