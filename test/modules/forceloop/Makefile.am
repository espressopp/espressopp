ESLIBS= \
        $(BOOST_PYTHON_LIBS) \
        $(top_builddir)/src/libespresso_common.la \
        $(top_builddir)/src/libgenlogger.la

check_PROGRAMS = forceloop benchmark

forceloop_SOURCES = main.cpp
forceloop_LDADD   = $(ESLIBS)

benchmark_SOURCES = benchmark.cpp
benchmark_LDADD   = $(ESLIBS)
