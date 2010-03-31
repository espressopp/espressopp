mpi4py_srcdir = mpi4py-1.2.1

noinst_HEADERS = \
	$(mpi4py_srcdir)/src/atimport.h \
	$(mpi4py_srcdir)/src/fallback.h \
	$(mpi4py_srcdir)/src/missing.h \
	$(mpi4py_srcdir)/src/vendor.h \
	$(mpi4py_srcdir)/src/dynload.h \
	$(mpi4py_srcdir)/src/config/unknown.h \
	$(mpi4py_srcdir)/src/config/mpich2.h \
	$(mpi4py_srcdir)/src/config/mpich2io.h \
	$(mpi4py_srcdir)/src/config/sgimpi.h \
	$(mpi4py_srcdir)/src/config/lammpi.h \
	$(mpi4py_srcdir)/src/config/lammpiio.h \
	$(mpi4py_srcdir)/src/config/mpich1.h \
	$(mpi4py_srcdir)/src/config/mpich1io.h \
	$(mpi4py_srcdir)/src/config/openmpi.h \
	$(mpi4py_srcdir)/src/compat/mpich1.h \
	$(mpi4py_srcdir)/src/compat/mpich2.h \
	$(mpi4py_srcdir)/src/compat/openmpi.h \
	$(mpi4py_srcdir)/src/compat/sgimpi.h \
	$(mpi4py_srcdir)/src/compat/sicortex.h \
	$(mpi4py_srcdir)/src/compat/lammpi.h \
	$(mpi4py_srcdir)/src/mpi4py.MPI.c \
	$(mpi4py_srcdir)/src/MPI.c \
	$(mpi4py_srcdir)/src/dynload.c

EXTRA_DIST = \
	$(mpi4py_srcdir)/LICENSE.txt \
	$(mpi4py_srcdir)/HISTORY.txt \
	$(mpi4py_srcdir)/README.txt \
	$(mpi4py_srcdir)/THANKS.txt

# Create a convenience library version of the python module,
# required for static linkage into pypresso
noinst_LTLIBRARIES = libmpi4py.la
libmpi4py_la_CPPFLAGS = -I $(srcdir)/$(mpi4py_srcdir)/src $(AM_CPPFLAGS)
libmpi4py_la_LDFLAGS = -no-undefined
libmpi4py_la_LIBADD = $(ESPRESSO_BOOST_LIBS)

libmpi4py_la_SOURCES = MPI.cpp dynload.cpp

# The python module MPI
# Build only if shared libs are built
if ENABLE_SHARED
python_LTLIBRARIES = MPI.la
MPI_la_LDFLAGS = -shared -module -no-undefined -avoid-version
MPI_la_LIBTOOLFLAGS = --tag=CXX
MPI_la_SOURCES =
MPI_la_LIBADD = libmpi4py.la $(ESPRESSO_BOOST_LIBS)
endif
