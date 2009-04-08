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
	estypes.hpp \
	Real3D.hpp \
	espresso_common.hpp

# the common part of all Espresso flavors
noinst_LTLIBRARIES = libespresso_common.la
libespresso_common_la_SOURCES= \
	Real3D.cpp \
	bindings.cpp \
	espresso_common.cpp
libespresso_common_la_LIBADD = \
	esutil/libesutil.la \
	hello/libhello.la \
	particles/libparticles.la \
	pairs/libpairs.la \
	interaction/libinteraction.la \
	bc/libbc.la \
	integrator/libintegrator.la \
	$(ESPRESSO_BOOST_LIBS) \
	$(PYTHON_LIBS) \
	$(MPI_LIBS) 

# The python module _espresso
# Build only if shared libs are built
if ENABLE_SHARED
python_LTLIBRARIES = _espresso.la
_espresso_la_SOURCES = _espresso.cpp
_espresso_la_LDFLAGS = -shared -module -avoid-version -no-undefined
_espresso_la_LIBADD = \
	libespresso_common.la \
	esutil/libpylogger.la
endif # ENABLE_SHARED

# the embedded python interpreter that makes _espresso and _boostmpi a builtin
bin_PROGRAMS = pypresso 
pypresso_SOURCES = pypresso.cpp
pypresso_LDADD = \
	libespresso_common.la \
	esutil/libpylogger.la \
	../contrib/boostmpi/libboostmpi.la \
	$(ESPRESSO_BOOST_LIBS) \
	$(PYTHON_LIBS) \
	$(MPI_LIBS)

##############################
## PYTHON
##############################

python_scripts=\
	__init__.py \
	eslogging.py \
	boostmpi.py \
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

