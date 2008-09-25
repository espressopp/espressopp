SUBDIRS = pmi espresso

lib_LTLIBRARIES = _hello.la

# TODO: results in a warning. Requires libltdl!

noinst_HEADERS = \
	hello.hpp \
	hello/HelloWorld.hpp

_hello_la_SOURCES = \
	hello/python_bindings.cpp \
	hello/HelloWorld.cpp 

# set the flags where to find PMI and BOOST
_hello_la_CPPFLAGS = $(BOOST_CPPFLAGS) 
_hello_la_LDFLAGS = -module $(BOOST_LDFLAGS)
_hello_la_LIBADD = pmi/libpmi.la $(BOOST_MPI_LIBS) $(BOOST_PYTHON_LIBS)
