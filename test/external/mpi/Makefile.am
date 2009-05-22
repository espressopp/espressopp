#
#  Makefile.am for automake
#
#  MPI test programs verify that MPI works correctly
#

noinst_PROGRAMS = probe barrier
check_PROGRAMS = cpi comm mpiTest

LDADD = $(MPI_LIBS)

cpi_SOURCES = cpi.cpp
comm_SOURCES = comm.cpp
mpiTest_SOURCES = mpiTest.cpp
probe_SOURCES = probe.cpp
probe_LDADD = $(LDADD) $(top_builddir)/src/esutil/Timer.lo
barrier_SOURCES = barrier.cpp
