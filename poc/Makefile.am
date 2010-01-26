python_SCRIPTS =
noinst_LTLIBRARIES =
noinst_PROGRAMS =
noinst_HEADERS =

include boost/SubMakefile.am
include boost.python/inheritance/SubMakefile.am
include boost.python/shared_ptr/SubMakefile.am
#include boost.python/keywords/SubMakefile.am
include boost.python/diamond/SubMakefile.am
include boost.python/real3D/SubMakefile.am
include forward_declaration/SubMakefile.am

include $(top_srcdir)/build-aux/macros/axes_python.am

noinst_PROGRAMS += wrappingReal3D
wrappingReal3D_SOURCES = wrappingReal3D.cpp
wrappingReal3D_LDADD = ../src/libespresso_common.la \
	$(ESPRESSO_BOOST_LIBS) $(PYTHON_LIBS)

noinst_PROGRAMS += inline
inline_SOURCES = inline.cpp
inline_LDADD = ../src/libespresso_common.la \
	$(ESPRESSO_BOOST_LIBS) $(PYTHON_LIBS)

noinst_PROGRAMS += packed_mpi
packed_mpi_SOURCES = packed_mpi.cpp
packed_mpi_LDADD = ../src/libespresso_common.la \
	$(ESPRESSO_BOOST_LIBS) $(PYTHON_LIBS)

noinst_PROGRAMS += boost_function
boost_function_SOURCES = boost_function.cpp
boost_function_LDADD = ../src/libespresso_common.la \
	$(ESPRESSO_BOOST_LIBS) $(PYTHON_LIBS)

noinst_PROGRAMS += interruptable_loop
interruptable_loop_SOURCES = interruptable_loop.cpp
interruptable_loop_LDADD = ../src/libespresso_common.la

noinst_PROGRAMS += crtp 
crtp_SOURCES = crtp.cpp

python_SCRIPTS += __init__.py
