# declare the subdirs
SUBDIRS = \
	include \
	esutil \
	hello \
	particles \
	pairs \
	interaction \
	bc \
	integrator
CLEANFILES =

# declare the headers
noinst_HEADERS = \
	espresso_common.hpp

# the common part of all Espresso flavors
noinst_LTLIBRARIES = libespresso_common.la
libespresso_common_la_SOURCES= \
	espresso_common.cpp
libespresso_common_la_LIBADD = \
	esutil/libesutil.la \
	hello/libhello.la \
	particles/libparticles.la \
	pairs/libpairs.la \
	interaction/libinteraction.la \
	bc/libbc.la \
	integrator/libintegrator.la \
	$(BOOST_MPI_LIBS) \
	$(BOOST_PYTHON_LIBS) \
	$(LOG4ESPP_LIBS)

# Rules to build the python module _espresso
python_LTLIBRARIES = _espresso.la
_espresso_la_SOURCES = _espresso.cpp
_espresso_la_LDFLAGS = -module -avoid-version -no-undefined
_espresso_la_LIBADD = libespresso_common.la

# the embedded python interpreter with static espresso library
bin_PROGRAMS = pypresso 
pypresso_SOURCES = pypresso.cpp
pypresso_LDADD = libespresso_common.la

##############################
## PYTHON
##############################

python_scripts=\
	__init__.py \
	esmpi.py \
	pmi.py \
	esutil.py \
	hello.py \
	particles.py \
	interaction.py

python_SCRIPTS=$(python_scripts)

# python script handling
include $(top_srcdir)/build-aux/macros/axes_python.am

clean-local: python-clean
all-local: python-all

