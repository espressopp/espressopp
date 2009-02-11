LIBS= \
	$(top_builddir)/src/cpp/particlestorage/libparticlestorage.la \
	$(top_builddir)/src/cpp/interaction/libinteraction.la \
	$(top_builddir)/src/cpp/pairs/libpairs.la \
	$(top_builddir)/src/cpp/esutil/libesutil.la \
	$(LOG4ESPP_LIBS)

check_PROGRAMS = forceloop benchmark

forceloop_SOURCES = main.cpp

forceloop_LDADD   = $(LIBS)

benchmark_SOURCES = benchmark.cpp
benchmark_LDADD   = $(LIBS)
