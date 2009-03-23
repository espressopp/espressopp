python_scripts=main.py hello.py

EXTRA_DIST = $(python_scripts)

hello_sources = HelloWorld.cpp HelloWorld.hpp

# build a (normally shared) library that will not be installed
check_LTLIBRARIES = _espresso.la

_espresso_la_SOURCES = _espresso.cpp $(hello_sources)
# the -module flag says that the library is to be used via dlopen
# -rpath is required by automake so that a shared library is actually built
_espresso_la_LDFLAGS = -module -avoid-version -no-undefined -rpath /nowhere
_espresso_la_LIBADD= $(BOOST_PYTHON_LIBS) $(PYTHON_LIBS)

# and a regular binary with embedded python
check_PROGRAMS = espresso

espresso_SOURCES= main_embedded.cpp _espresso.cpp $(hello_sources)
# work around a bug in automake that does not allow to use the same
# object in a library and a binary. We also cannot link against
# _espresso.so, because that causes the shared library to be linked
# against
espresso_CPPFLAGS = $(AM_CPPFLAGS)
espresso_LDADD=$(BOOST_PYTHON_LIBS) $(PYTHON_LIBS)

# python script handling
include $(top_srcdir)/build-aux/macros/axes_python.am

all-local: python-all
clean-local: python-clean
