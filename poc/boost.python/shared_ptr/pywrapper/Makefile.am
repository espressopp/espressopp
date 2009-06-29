SIMULATION=$(top_srcdir)/shared_ptr/cpp
SIMULATIONLIB=$(top_builddir)/shared_ptr/cpp/libsimulation.la

AM_CPPFLAGS+=-I$(SIMULATION)

LIBS=$(BOOST_PYTHON_LIBS)

python_scripts = cSimulation.py Simulation.py wrapper.py

# build a convenience library that will not be installed

lib_LTLIBRARIES = _Simulation.la

# Note: flag -module is needed to allow for a library without prefix lib

_Simulation_la_LDFLAGS = -module
_Simulation_la_SOURCES = wrapper.cpp
_Simulation_la_LIBADD  = $(SIMULATIONLIB) 

include $(top_srcdir)/../config/axes_dllfixup.am
include $(top_srcdir)/../config/axes_python.am

all-local: dllfixup python-all 

clean-local: dllclean python-clean
