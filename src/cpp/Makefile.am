# declare the subdirs
SUBDIRS = logging hello
EXTRA_DIST =
CLEANFILES =

if HAVE_MPI
  SUBDIRS += mpi pmi
else
  EXTRA_DIST += mpi pmi
endif

# declare the compile flags
AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(LOG4ESPP_CPPFLAGS)
AM_LDFLAGS = $(BOOST_LDFLAGS)

# declare the headers
noinst_HEADERS = \
	espresso_common.hpp \
	logging.hpp \
	hello.hpp

if HAVE_MPI
noinst_HEADERS += mpi.hpp
endif

# the common part of all Espresso flavors
noinst_LTLIBRARIES = libespresso_common.la
libespresso_common_la_SOURCES= espresso_common.cpp
libespresso_common_la_LIBADD = \
	hello/libhello.la \
	logging/liblogging.la \
	interaction/libinteraction.la \
	$(BOOST_PYTHON_LIBS) $(LOG4ESPP_LIBS)

if HAVE_MPI
libespresso_common_la_LIBADD += \
	pmi/libpmi.la \
	mpi/libmpi.la \
	$(BOOST_MPI_LIBS)
endif

# the program wrappers
bin_PROGRAMS = espresso # espresso_worker

espresso_SOURCES = ../cpp/espresso.cpp
espresso_LDADD = ../cpp/libespresso_common.la

# espresso_worker_SOURCES = ../cpp/espresso_main.cpp
# espresso_worker_CPPFLAGS = -DWORKER_ONLY
# espresso_worker_LDADD = ../cpp/libespresso_common.la

# Rules to build the python module _espresso
python_LTLIBRARIES = _espresso.la
_espresso_la_SOURCES = _espresso.cpp
_espresso_la_LDFLAGS = -module -no-undefined 
_espresso_la_LIBADD = libespresso_common.la

