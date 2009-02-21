LIBS= \
	$(top_builddir)/src/particlestorage/libparticlestorage.la \
	$(top_builddir)/src/interaction/libinteraction.la \
        $(top_builddir)/src/bc/libbc.la \
	$(top_builddir)/src/pairs/libpairs.la \
	$(top_builddir)/src/integrator/libintegrator.la \
	$(top_builddir)/src/esutil/libesutil.la \
        $(BOOST_PYTHON_LIBS) \
        $(top_builddir)/src/libespresso_common.la
	$(LOG4ESPP_LIBS)

if HAVE_MPI
LIBS +=\
	$(top_builddir)/src/pmi/libpmi.la \
        $(BOOST_MPI_LIBS)
endif

check_PROGRAMS = forceloop benchmark

forceloop_SOURCES = main.cpp
forceloop_LDADD   = $(LIBS)

benchmark_SOURCES = benchmark.cpp
benchmark_LDADD   = $(LIBS)
