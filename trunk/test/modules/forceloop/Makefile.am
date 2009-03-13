LIBS= \
        $(BOOST_PYTHON_LIBS) \
        $(top_builddir)/src/libespresso_common.la
	$(LOG4ESPP_LIBS)

check_PROGRAMS = forceloop benchmark

forceloop_SOURCES = main.cpp
forceloop_LDADD   = $(LIBS)

benchmark_SOURCES = benchmark.cpp
benchmark_LDADD   = $(LIBS)
