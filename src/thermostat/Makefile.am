subpkgincludedir=$(pkgincludedir)/thermostat

subpkginclude_HEADERS = \
        bindings.hpp \
        Thermostat.hpp \
	Langevin.hpp 

noinst_LTLIBRARIES = libthermostat.la

libthermostat_la_SOURCES = \
        bindings.cpp \
        Thermostat.cpp \
	Langevin.cpp

include $(top_srcdir)/build-aux/macros/axes_python.am
#python_SCRIPTS = __init__.py
