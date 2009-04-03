DIST_SUBDIRS=unittest

subpkgincludedir=$(pkgincludedir)/particles

subpkginclude_HEADERS = \
	bindings.hpp \
	Computer.hpp \
	PythonComputer.hpp \
	ParticleHandle.hpp \
	Property.hpp \
	Storage.hpp \
	Set.hpp \
	All.hpp

noinst_LTLIBRARIES = libparticles.la

libparticles_la_SOURCES = \
	PythonComputer.cpp \
	bindings.cpp \
	Storage.cpp \
	All.cpp
