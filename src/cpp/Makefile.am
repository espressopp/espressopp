SUBDIRS = pmi espresso

# make the dll's in this directory available even in the
# build tree
include $(top_srcdir)/build-aux/macros/axes_dllfixup.am
all-local: dllfixup espresso.py
#FIXME: espresso.py and the rest of python needs to be installed also
clean-local: dllclean

AM_CPPFLAGS = $(BOOST_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)
LIBS = $(BOOST_PYTHON_LIBS) $(BOOST_MPI_LIBS)
LIBADD = pmi/libpmi.la $(BOOST_MPI_LIBS) $(BOOST_PYTHON_LIBS)

lib_LTLIBRARIES = _hello.la

# TODO: results in a warning. Requires libltdl!

noinst_HEADERS = \
	hello.hpp \
	hello/HelloWorld.hpp

_hello_la_SOURCES = \
	hello/python_bindings.cpp \
	hello/HelloWorld.cpp 

_hello_la_LDFLAGS = -module -no-undefined 
