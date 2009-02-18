# declare the subdirs
SUBDIRS = include esutil
EXTRA_DIST =
CLEANFILES =

if HAVE_MPI
  SUBDIRS += pmi
else
  EXTRA_DIST += pmi
endif

SUBDIRS += \
	hello interaction integrator \
	pairs bc particleset particlestorage


# declare the headers
noinst_HEADERS = \
	espresso_common.hpp

# the common part of all Espresso flavors
noinst_LTLIBRARIES = libespresso_common.la
libespresso_common_la_SOURCES= espresso_common.cpp
libespresso_common_la_LIBADD = \
	hello/libhello.la \
	interaction/libinteraction.la \
	$(BOOST_PYTHON_LIBS) $(LOG4ESPP_LIBS)

if HAVE_MPI
libespresso_common_la_LIBADD += \
	pmi/libpmi.la \
	$(BOOST_MPI_LIBS)
endif

# the program wrappers
bin_PROGRAMS = espresso # espresso_worker

espresso_SOURCES = espresso.cpp
espresso_LDADD = libespresso_common.la

# espresso_worker_SOURCES = ../cpp/espresso_main.cpp
# espresso_worker_CPPFLAGS = -DWORKER_ONLY
# espresso_worker_LDADD = ../cpp/libespresso_common.la

# Rules to build the python module _espresso
python_LTLIBRARIES = _espresso.la
_espresso_la_SOURCES = _espresso.cpp
_espresso_la_LDFLAGS = -module -no-undefined 
_espresso_la_LIBADD = libespresso_common.la

python_scripts=hello.py 

## THE SCRIPTS
python_SCRIPTS=$(python_scripts)

# python script handling
include $(top_srcdir)/build-aux/macros/axes_python.am

# clean the __init__.py-files from the build dir
MOSTLYCLEANFILES=__init__.py __init__.pyc __init__.install.py

# rename __init__.build.py to __init__.py in the build dir
# copy __init__.install.py into the build dir
all-local: python-all
	sed -e "s,[@]srcdir[@],'$(abs_top_builddir)/src',g" \
	$(srcdir)/__init__.build.py > __init__.py

# install __init__.install.py into pkglibdir, renaming it to __init__.py
install-exec-local:
	$(install_sh_SCRIPT) $(srcdir)/__init__.install.py $(DESTDIR)$(pythondir)/__init__.py

# remove __init__.py and __init__.pyc from pkglib
uninstall-local:
	-rm $(DESTDIR)$(pythondir)/__init__.py
	-rm $(DESTDIR)$(pythondir)/__init__.pyc


