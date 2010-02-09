nobase_pkginclude_HEADERS =
python_SCRIPTS =
noinst_LTLIBRARIES =
libespresso_common_la_SOURCES =
DIST_SUBDIRS = 
EXTRA_DIST =
python_LTLIBRARIES =

include SubMakefile.am
include include/SubMakefile.am
include esutil/SubMakefile.am
include bc/SubMakefile.am
include storage/SubMakefile.am
include iterator/SubMakefile.am
include interaction/SubMakefile.am
include integrator/SubMakefile.am

CLEANFILES =

# declare the headers
noinst_HEADERS = \
	espresso_common.hpp

# the common part of all Espresso flavors
noinst_LTLIBRARIES += libespresso_common.la
libespresso_common_la_SOURCES += \
	espresso_common.cpp
libespresso_common_la_LIBADD = \
	$(ESPRESSO_BOOST_LIBS) 

# The python module _espresso
# Build only if shared libs are built
if ENABLE_SHARED
python_LTLIBRARIES += _espresso.la
_espresso_la_SOURCES = _espresso.cpp
_espresso_la_LDFLAGS = -shared -module -avoid-version -no-undefined
_espresso_la_LIBADD = \
	libespresso_common.la \
	libpylogger.la
endif # ENABLE_SHARED

# the embedded python interpreter that makes _espresso and _boostmpi a builtin
bin_PROGRAMS = pypresso 
pypresso_SOURCES = pypresso.cpp
pypresso_LDADD = \
	libespresso_common.la \
	../contrib/mpi4py/libmpi4py.la \
	libpylogger.la \
	$(PYTHON_LIBS) \
	$(MPI_LIBS)

# check whether a link to pypresso in the bin directory exists
all-local:
	if test -x $(top_builddir)/src/pypresso && \
	    test ! -L $(top_builddir)/bin/pypresso; then \
	  $(LN_S) $(top_builddir)/src/pypresso $(top_builddir)/bin/; \
	fi

##############################
## PYTHON
##############################

# python script handling
include $(top_srcdir)/build-aux/macros/axes_python.am

install-exec-local:
# install espresso setup scripts
	$(install_sh_SCRIPT) __setupPath.py.install $(DESTDIR)$(pkglibdir)/__setupPath.py

uninstall-local:
	rm -f $(DESTDIR)$(pkglibdir)/__setupPath.py

