#
#  Makefile.am for automatake
#
#  Test programs for the kernel
#

LIBS       = $(LOG4ESPP_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

check_PROGRAMS = BasicPropertyTest ParticleContainerTest

TESTS = $(check_PROGRAMS)

BasicPropertyTest_SOURCES = BasicPropertyTest.cpp ParticleContainerMock.cpp
BasicPropertyTest_CPPFLAGS = $(BOOST_CPPFLAGS) $(LOG4ESPP_CPPFLAGS)
BasicPropertyTest_LDFLAGS = $(BOOST_LDFLAGS) $(LOG4ESPP_LDFLAGS)
BasicPropertyTest_LDADD = ../BasicProperty.o

ParticleContainerTest_SOURCES  = ParticleContainerTest.cpp 
ParticleContainerTest_CPPFLAGS = $(BOOST_CPPFLAGS) $(LOG4ESPP_CPPFLAGS)
ParticleContainerTest_LDFLAGS  = $(BOOST_LDFLAGS) $(LOG4ESPP_LDFLAGS)
ParticleContainerTest_LDADD    = ../BasicProperty.o ../ParticleContainer.o

INCLUDES = -I../
