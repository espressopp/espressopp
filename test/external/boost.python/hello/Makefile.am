AM_CPPFLAGS=$(BOOST_CPPFLAGS)
AM_LDFLAGS=$(BOOST_LDFLAGS)
LIBS= $(BOOST_PYTHON_LIBS)

libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# build a convenience library that will not be installed
noinst_LTLIBRARIES = libhello.la

libhello_la_LDFLAGS = -no-undefined
libhello_la_SOURCES = python_bindings.cpp \
	HelloWorld.hpp HelloWorld.cpp