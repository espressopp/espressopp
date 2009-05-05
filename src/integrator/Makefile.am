subpkgincludedir=$(pkgincludedir)/integrator

subpkginclude_HEADERS = \
        bindings.hpp\
	MDIntegrator.hpp \
	VelocityVerlet.hpp 

noinst_LTLIBRARIES = libintegrator.la

libintegrator_la_SOURCES = \
        bindings.cpp\
	MDIntegrator.cpp \
	VelocityVerlet.cpp

include $(top_srcdir)/build-aux/macros/axes_python.am
python_SCRIPTS = __init__.py
