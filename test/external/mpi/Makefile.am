#
#  Makefile.am for automake
#
#  MPI test programs verify that MPI works correctly
#

check_PROGRAMS = cpi comm probe mpiTest

LDADD = $(MPI_LIBS)

cpi_SOURCES = cpi.cpp
comm_SOURCES = comm.cpp
mpiTest_SOURCES = mpiTest.cpp
probe_SOURCES = probe.cpp
probe_LDADD = $(MPI_LIBS) $(top_builddir)/src/esutil/libesutil.la

