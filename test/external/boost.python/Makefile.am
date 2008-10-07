AM_CPPFLAGS=$(BOOST_CPPFLAGS)
AM_LDFLAGS=$(BOOST_LDFLAGS)
LIBS= $(BOOST_PYTHON_LIBS)

python_scripts=main.py hello.py

EXTRA_DIST = $(python_scripts)

# build a (normally shared) library that will be installed
check_LTLIBRARIES = _espresso.la

# and a regular binary with embedded python
check_PROGRAMS = espresso

hello_sources = HelloWorld.cpp HelloWorld.hpp

# the -module flag says that the library is to be used via dlopen
_espresso_la_LDFLAGS = -module -no-undefined
_espresso_la_SOURCES = _espresso.cpp $(hello_sources)

espresso_SOURCES= main_embedded.cpp _espresso.cpp $(hello_sources)
# work around a bug in automake that does not allow to use the same
# object in a library and a binary. We also cannot link against
# _espresso.so, because that causes the shared library to be linked
# against
espresso_CPPFLAGS = $(AM_CPPFLAGS)

# python script handling
include $(top_srcdir)/build-aux/macros/axes_python.am

# make the dll's in this directory available even in the
# build tree
include $(top_srcdir)/build-aux/macros/axes_dllfixup.am

all-local: dllfixup python-all
clean-local: dllclean python-clean
