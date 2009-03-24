# declare the subdirs
SUBDIRS = \
	include \
	esutil \
	base \
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
	base/libbase.la \
	hello/libhello.la \
	particles/libparticles.la \
	pairs/libpairs.la \
	interaction/libinteraction.la \
	bc/libbc.la \
	integrator/libintegrator.la \
	$(PYTHON_LIBS) \
	$(MPI_LIBS) \
	$(BOOST_SERIALIZATION_LIBS) \
	$(BOOST_MPI_LIBS) \
	$(BOOST_PYTHON_LIBS) \
	$(BOOST_MPI_PYTHON_LIBS) \
	$(LOG4ESPP_LIBS)

# The python module _espresso
# Build only if shared libs are built
if ENABLE_SHARED
python_LTLIBRARIES = _espresso.la
_espresso_la_SOURCES = _espresso.cpp
_espresso_la_LDFLAGS = -shared -module -avoid-version -no-undefined
_espresso_la_LIBADD = libespresso_common.la
endif # ENABLE_SHARED

# the embedded python interpreter with static espresso library
bin_PROGRAMS = pypresso 
pypresso_SOURCES = pypresso.cpp
pypresso_LDADD = libespresso_common.la \
	$(PYTHON_LIBS) \
	$(MPI_LIBS) \
	$(BOOST_SERIALIZATION_LIBS) \
	$(BOOST_MPI_LIBS) \
	$(BOOST_PYTHON_LIBS) \
	$(BOOST_MPI_PYTHON_LIBS) \
	$(LOG4ESPP_LIBS)


##############################
## PYTHON
##############################

python_scripts=\
	__init__.py \
	base.py \
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

