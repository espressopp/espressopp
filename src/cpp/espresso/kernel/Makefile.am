SUBDIRS = test

## Open problem: enabling of gcov at a higher level

# AM_CPPFLAGS = -fprofile-arcs -ftest-coverage
# AM_LDFLAGS  = -fprofile-arcs -ftest-coverage

lib_LTLIBRARIES = libesppkernel.la

libesppkernel_la_SOURCES = BasicProperty.cpp BasicProperty.hpp \
                           ParticleContainer.cpp ParticleContainer.hpp 

#  for clean: remove all the files generate by gcov code coverage
#

CLEANFILES = *.gcov *.gcda *.gcno
