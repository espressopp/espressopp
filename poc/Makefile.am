python_SCRIPTS =
noinst_LTLIBRARIES =
noinst_PROGRAMS =

include boost/SubMakefile.am
include boost.python/inheritance/SubMakefile.am
include boost.python/shared_ptr/SubMakefile.am
#include boost.python/keywords/SubMakefile.am
include boost.python/diamond/SubMakefile.am

include $(top_srcdir)/build-aux/macros/axes_python.am

noinst_PROGRAMS += inline
inline_SOURCES = inline.cpp
inline_LDADD = ../src/libespresso_common.la \
	$(ESPRESSO_BOOST_LIBS) $(PYTHON_LIBS)

noinst_PROGRAMS += boost_function
boost_function_SOURCES = boost_function.cpp
boost_function_LDADD = ../src/libespresso_common.la \
	$(ESPRESSO_BOOST_LIBS) $(PYTHON_LIBS)

noinst_PROGRAMS += crtp
crtp_SOURCES = crtp.cpp

python_SCRIPTS += __init__.py
