DIST_SUBDIRS=unittest

subpkgincludedir=$(pkgincludedir)/particles

subpkginclude_HEADERS = \
	bindings.hpp \
	Computer.hpp \
	PythonComputer.hpp \
	ParticleHandle.hpp \
	PropertyHandle.hpp \
	Storage.hpp \
	Set.hpp \
	All.hpp

noinst_LTLIBRARIES = libparticles.la

libparticles_la_SOURCES = \
	PythonComputer.cpp \
	bindings.cpp \
	Storage.cpp \
	All.cpp

include $(top_srcdir)/build-aux/macros/axes_python.am
python_SCRIPTS = __init__.py
