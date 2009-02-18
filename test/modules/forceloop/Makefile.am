LIBS= \
	$(top_builddir)/src/particlestorage/libparticlestorage.la \
	$(top_builddir)/src/interaction/libinteraction.la \
	$(top_builddir)/src/pairs/libpairs.la \
	$(top_builddir)/src/integrator/libintegrator.la \
	$(top_builddir)/src/pmi/libpmi.la \
	$(top_builddir)/src/esutil/libesutil.la \
        $(BOOST_MPI_LIBS) \
        $(BOOST_PYTHON_LIBS) \
        $(top_builddir)/src/libespresso_common.la
	$(LOG4ESPP_LIBS)

check_PROGRAMS = forceloop benchmark

forceloop_SOURCES = main.cpp

forceloop_LDADD   = $(LIBS)

benchmark_SOURCES = benchmark.cpp
benchmark_LDADD   = $(LIBS)
