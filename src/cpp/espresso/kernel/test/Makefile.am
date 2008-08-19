#
#  Makefile.am for automatake
#
#  Test programs for the kernel
#

check_PROGRAMS = BasicPropertyTest ParticleContainerTest

TESTS = $(check_PROGRAMS)

BasicPropertyTest_SOURCES = BasicPropertyTest.cpp ParticleContainerMock.cpp
BasicPropertyTest_CPPFLAGS = $(BOOST_CPPFLAGS)
BasicPropertyTest_LDFLAGS = $(BOOST_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)
BasicPropertyTest_LDADD = ../BasicProperty.o

ParticleContainerTest_SOURCES  = ParticleContainerTest.cpp 
ParticleContainerTest_CPPFLAGS = $(BOOST_CPPFLAGS)
ParticleContainerTest_LDFLAGS  = $(BOOST_LDFLAGS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)
ParticleContainerTest_LDADD    = ../BasicProperty.o ../ParticleContainer.o

INCLUDES = -I$(top_srcdir)/src/cpp/kernel
