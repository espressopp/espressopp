AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(LOG4ESPP_CPPFLAGS)
AM_LDFLAGS  = $(BOOST_LDFLAGS)

PARTICLESTORAGE_LIB=$(top_builddir)/src/cpp/particlestorage/libparticlestorage.la
INTERACTION_LIB=$(top_builddir)/src/cpp/interaction/libinteraction.la

bin_PROGRAMS = forceloop

forceloop_SOURCES = main.cpp

forceloop_LDADD   = $(PARTICLESTORAGE_LIB) $(INTERACTION_LIB) $(LOG4ESPP_LIBS)
