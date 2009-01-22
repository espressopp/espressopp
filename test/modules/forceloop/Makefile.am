AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(LOG4ESPP_CPPFLAGS)
AM_LDFLAGS  = $(BOOST_LDFLAGS)

PARTICLESTORAGE_LIB=$(top_builddir)/src/cpp/particlestorage/libparticlestorage.la

bin_PROGRAMS = forceloop

forceloop_SOURCES = main.cpp

forceloop_LDADD   = $(PARTICLESTORAGE_LIB)
