LIBS= \
	$(top_builddir)/src/particlestorage/libparticlestorage.la \
	$(top_builddir)/src/interaction/libinteraction.la \
	$(top_builddir)/src/pairs/libpairs.la \
	$(top_builddir)/src/integrator/libintegrator.la \
	$(top_builddir)/src/esutil/libesutil.la \
	$(LOG4ESPP_LIBS)

check_PROGRAMS = forceloop benchmark

forceloop_SOURCES = main.cpp

forceloop_LDADD   = $(LIBS)

benchmark_SOURCES = benchmark.cpp
benchmark_LDADD   = $(LIBS)
