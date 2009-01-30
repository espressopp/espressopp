PARTICLESTORAGE_LIB=$(top_builddir)/src/cpp/particlestorage/libparticlestorage.la
INTERACTION_LIB=$(top_builddir)/src/cpp/interaction/libinteraction.la
PAIRS_LIB=$(top_builddir)/src/cpp/pairs/libpairs.la

bin_PROGRAMS = forceloop

forceloop_SOURCES = main.cpp

forceloop_LDADD   = $(PARTICLESTORAGE_LIB) $(INTERACTION_LIB) $(PAIRS_LIB) $(LOG4ESPP_LIBS)
