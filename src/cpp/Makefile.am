SUBDIRS = espresso

lib_LTLIBRARIES = _hello.la

# TODO: results in a warning. Requires libltdl!

noinst_HEADERS = \
	hello.hpp \
	hello/HelloWorld.hpp

_hello_la_SOURCES = \
	hello/python_bindings.cpp \
	hello/HelloWorld.cpp 

# set the flags where to find PMI and BOOST
_hello_la_CPPFLAGS = $(BOOST_CPPFLAGS) -I$(top_srcdir)/contrib/pmi
#_hello_la_LIBADD = $(top_builddir)/contrib/pmi/pmi/libpmi.la
_hello_la_LDFLAGS = -module
