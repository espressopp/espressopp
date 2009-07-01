noinst_HEADERS = \
	include/boostmpi/config.hpp \
	include/boostmpi/serialize.hpp \
	include/boostmpi/skeleton_and_content.hpp \
	include/boostmpi.hpp \
	src/request_with_value.hpp

# Create a convenience library version of the python module, 
# required for static linkage into pypresso
noinst_LTLIBRARIES = libboostmpi.la
libboostmpi_la_LIBTOOLFLAGS = --tag=CXX
libboostmpi_la_LDFLAGS = -no-undefined
libboostmpi_la_LIBADD = \
	$(ESPRESSO_BOOST_LIBS) \
	$(MPI_LIBS) $(PYTHON_LIBS) 

libboostmpi_la_SOURCES = \
	src/collectives.cpp \
	src/datatypes.cpp \
	src/documentation.cpp \
	src/module.cpp \
	src/py_communicator.cpp \
	src/py_environment.cpp \
	src/py_exception.cpp \
	src/py_nonblocking.cpp \
	src/py_request.cpp \
	src/py_timer.cpp \
	src/serialize.cpp \
	src/skeleton_and_content.cpp \
	src/status.cpp

# The python module boostmpi
# Build only if shared libs are built
if ENABLE_SHARED
python_LTLIBRARIES = _boostmpi.la
_boostmpi_la_LDFLAGS = -shared -module -no-undefined -avoid-version
_boostmpi_la_SOURCES =
_boostmpi_la_LIBADD = \
	libboostmpi.la
_boostmpi_la_LIBS = \
	$(ESPRESSO_BOOST_LIBS) \
	$(MPI_LIBS) $(PYTHON_LIBS) 
endif

