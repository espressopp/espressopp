#
#  Makefile.am for automake
#
#  C++ main program using the Espresso++ kernel
#

bin_PROGRAMS = container halley

INCLUDES = -I$(top_srcdir)/src/cpp/espresso/kernel 

LDADD = $(top_builddir)/src/cpp/espresso/kernel/libesppkernel.la

AM_INCLUDE = $(LOG4ESPP_CPPFLAGS)
AM_LDFLAGS = $(LOG4ESPP_LDFLAGS)

LIBS = $(LOG4ESPP_LIBS)

container_SOURCES = container.cpp

halley_SOURCES = halley.cpp

EXTRA_DIST=log4cxx.properties log4cpp.init

