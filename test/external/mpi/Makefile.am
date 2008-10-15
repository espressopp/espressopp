#
#  Makefile.am for automake
#
#  MPI test programs verify that MPI works correctly
#

check_PROGRAMS = cpi comm mpiTest

LIBS = $(MPI_LIBS)

cpi_SOURCES = cpi.cpp
comm_SOURCES = comm.cpp
mpiTest_SOURCES = mpiTest.cpp

check-local:
	$(top_builddir)/bin/espp_mpiexec.sh -n 2 mpiTest
	$(top_builddir)/bin/espp_mpiexec.sh -n 4 mpiTest
