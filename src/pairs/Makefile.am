subpkgincludedir=$(pkgincludedir)/pairs

subpkginclude_HEADERS = \
	All.hpp \
	PairEnergyComputer.hpp \
	ParticlePairComputer.hpp \
	PairComputerList.hpp \
	PairForceComputer.hpp \
	ParticlePairs.hpp

noinst_LTLIBRARIES = libpairs.la

libpairs_la_SOURCES = \
	All.cpp \
        List.cpp


check_PROGRAMS = ListTest 

TESTS = $(check_PROGRAMS)

ListTest_SOURCES = unittest/ListTest.cpp List.cpp \
                   ../particlestorage/ParticleStorage.cpp \
                   ../esutil/TupleVector.cpp \
                   ../bc/PBC.cpp 

ListTest_CPPFLAGS = $(BOOST_CPPFLAGS)
ListTest_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

