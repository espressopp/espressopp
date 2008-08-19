SUBDIRS = espresso

lib_LTLIBRARIES = _hello.la

# TODO: results in a warning. Requires libltdl!
_hello_la_SOURCES = \
	hello.hpp \
	hello/python_bindings.cpp \
	hello/HelloController.cpp hello/HelloController.hpp
_hello_la_LDFLAGS = -module
