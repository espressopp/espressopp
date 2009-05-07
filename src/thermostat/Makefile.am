subpkgincludedir=$(pkgincludedir)/thermostat

subpkginclude_HEADERS = \
        bindings.hpp\
	Langevin.hpp 

noinst_LTLIBRARIES = libthermostat.la

libthermostat_la_SOURCES = \
        bindings.cpp\
	Langevin.cpp

include $(top_srcdir)/build-aux/macros/axes_python.am
#python_SCRIPTS = __init__.py
