SUBDIRS = test

## Open problem: enabling of gcov at a higher level

# AM_CPPFLAGS = -fprofile-arcs -ftest-coverage
# AM_LDFLAGS  = -fprofile-arcs -ftest-coverage

lib_LTLIBRARIES = libesppkernel.la

# we compile everything with mpicxx so this is not necesssary, but just in case
# to remember how it could be done with cxx 
# LIBS = -L$(MPI_HOME)/lib -lmpich -lpthread -lrt
# INCLUDES = -I$(MPI_HOME)/include

libesppkernel_la_SOURCES = BasicProperty.cpp BasicProperty.hpp \
                           ParticleContainer.cpp ParticleContainer.hpp 

#  for clean: remove all the files generate by gcov code coverage
#

CLEANFILES = *.gcov *.gcda *.gcno
