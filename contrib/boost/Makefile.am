noinst_LTLIBRARIES=
python_LTLIBRARIES=

# This library bundles all boost libraries needed by espresso
pkglib_LTLIBRARIES = libespresso_boost.la
libespresso_boost_la_LDFLAGS = -no-undefined
libespresso_boost_la_SOURCES =
libespresso_boost_la_LIBADD =

# These makefiles add sources and libs to libespresso_boost.la
include Makefile-headers.am
include Makefile-serialization.am
include Makefile-mpi.am
include Makefile-python.am

# This builds the unit test framework
include Makefile-test.am

