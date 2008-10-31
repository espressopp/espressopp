#
#  Makefile.am for automake
#
#  MPI test programs verify that MPI works correctly
#

noinst_PROGRAMS = cpi comm mpiTest

LIBS = $(MPI_LIBS)

cpi_SOURCES = cpi.cpp
comm_SOURCES = comm.cpp
mpiTest_SOURCES = mpiTest.cpp

.PHONY: check-extended

localtest= $(top_builddir)/bin/runtests.sh -tests $(top_srcdir)/test/tests.txt -fromdir $(srcdir) 
check: $(ESPP_TESTS)
	$(localtest)
check-extended: $(ESPP_TESTS)
	$(localtest) -notag extended
