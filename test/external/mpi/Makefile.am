#
#  Makefile.am for automake
#
#  MPI test programs verify that MPI works correctly
#

noinst_PROGRAMS = cpi comm

check_PROGRAMS = mpiTest

AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(LOG4ESPP_CPPFLAGS)
AM_LDFLAGS  = $(BOOST_LDFLAGS) $(LOG4ESPP_LDFLAGS)

cpi_SOURCES = cpi.cpp
comm_SOURCES = comm.cpp

mpiTest_SOURCES = mpiTest.cpp

check-local:
	$(top_builddir)/bin/espp_mpiexec.sh -n 2 mpiTest
	$(top_builddir)/bin/espp_mpiexec.sh -n 4 mpiTest
