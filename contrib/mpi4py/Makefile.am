noinst_HEADERS = \
	atimport.h \
	fallback.h \
	missing.h \
	vendor.h \
	dynload.h \
	config/unknown.h \
	config/mpich2.h \
	config/mpich2io.h \
	config/sgimpi.h \
	config/lammpi.h \
	config/lammpiio.h \
	config/mpich1.h \
	config/mpich1io.h \
	config/openmpi.h \
	compat/mpich1.h \
	compat/mpich2.h \
	compat/openmpi.h \
	compat/sgimpi.h \
	compat/sicortex.h \
	compat/lammpi.h \
	mpi4py.MPI.c

# Create a convenience library version of the python module,
# required for static linkage into pypresso
CC=$(CXX)
noinst_LTLIBRARIES = libmpi4py.la
libmpi4py_la_LDFLAGS = -no-undefined
libmpi4py_la_LIBADD = $(ESPRESSO_BOOST_LIBS)

libmpi4py_la_SOURCES = MPI.c dynload.c

# The python module boostmpi
# Build only if shared libs are built
if ENABLE_SHARED
python_LTLIBRARIES = MPI.la
MPI_la_LDFLAGS = -shared -module -no-undefined -avoid-version
MPI_la_LIBTOOLFLAGS = --tag=CXX
MPI_la_SOURCES =
MPI_la_LIBADD = libmpi4py.la $(ESPRESSO_BOOST_LIBS)
endif
