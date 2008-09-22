#
#  Makefile.am for automake
#
#  C++ main program using the Espresso++ kernel
#

bin_PROGRAMS = container halley

INCLUDES = -I$(top_srcdir)/src/cpp/espresso/kernel 

LDADD = $(top_builddir)/src/cpp/espresso/kernel/libesppkernel.la

LIBS = $(LOG4ESPP_LIBS)

container_SOURCES = container.cpp
container_CPPFLAGS = $(LOG4ESPP_CPPFLAGS)
container_LDFLAGS = $(LOG4ESPP_LDFLAGS)

halley_SOURCES = halley.cpp
halley_CPPFLAGS = $(LOG4ESPP_CPPFLAGS)
halley_LDFLAGS = $(LOG4ESPP_LDFLAGS)

EXTRA_DIST=log4cxx.properties log4cpp.init

