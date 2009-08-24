LIBS= \
        $(BOOST_PYTHON_LIBS) \
        $(top_builddir)/src/libespresso_common.la \
        $(top_builddir)/src/libgenlogger.la

check_PROGRAMS = blockvector

blockvector_SOURCES = blockvector.cpp
blockvector_LDADD   = $(LIBS)
