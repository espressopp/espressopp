#
#  Makefile.am for automake
#
#  C++ main program using the Espresso++ kernel
#

bin_PROGRAMS = container

container_SOURCES = container.cpp
container_CPPFLAGS = $(LOG4ESPP_CPPFLAGS)
container_LDFLAGS = $(LOG4ESPP_LDFLAGS)

INCLUDES = -I$(top_srcdir)/src/cpp/espresso/kernel 

LDADD = $(top_srcdir)/src/cpp/espresso/kernel/libesppkernel.la

LIBS = $(LOG4ESPP_LIBS)
