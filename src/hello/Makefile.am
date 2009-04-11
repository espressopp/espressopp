subpkgincludedir=$(pkgincludedir)/hello
subpkginclude_HEADERS = \
	bindings.hpp \
	HelloWorld.hpp

noinst_LTLIBRARIES=libhello.la
libhello_la_SOURCES = \
	bindings.cpp \
	HelloWorld.cpp

include $(top_srcdir)/build-aux/macros/axes_python.am
python_SCRIPTS = __init__.py
