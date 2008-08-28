libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

# this is the final library, but packed as convenience library,
# since we want to version it the same way as the boost it was
# compiled against.
noinst_LTLIBRARIES = libboost_mpi.la
libboost_mpi_la_LDFLAGS = -no-undefined
libboost_mpi_la_SOURCES = \
    broadcast.cpp \
    communicator.cpp \
    computation_tree.cpp \
    content_oarchive.cpp \
    environment.cpp \
    exception.cpp \
    graph_communicator.cpp \
    group.cpp \
    intercommunicator.cpp \
    mpi_datatype_cache.cpp \
    mpi_datatype_oarchive.cpp \
    packed_iarchive.cpp \
    packed_oarchive.cpp \
    packed_skeleton_iarchive.cpp \
    packed_skeleton_oarchive.cpp \
    point_to_point.cpp \
    request.cpp \
    text_skeleton_oarchive.cpp

# Manually build the versioned library, since automake/libtool
# don't support determining the version only during configuration

VERSIONLIB=libboost_mpi@BOOST_SUFFIX@.la

$(VERSIONLIB): libboost_mpi.la
	$(LIBTOOL) --tag=CXX $(AM_LIBTOOLFLAGS) \
	        $(LIBTOOLFLAGS) --mode=link $(CXXLD) $(AM_CXXFLAGS) \
		$(CXXFLAGS) $(libboost_mpi_la_LDFLAGS) $(LDFLAGS) -o $@ \
		-rpath $(libdir) libboost_mpi.la $(LIBS)

install-versionlib: versionlib
	$(LIBTOOL) --mode=install $(INSTALL) $(INSTALL_STRIP_FLAG) $(VERSIONLIB) "$(DESTDIR)$(libdir)";


.PHONY: versionlib install-versionlib uninstall-versionlib clean-versionlib

all-local: $(VERSIONLIB)
install-exec-local: install-versionlib
uninstall-local: uninstall-versionlib
clean-local: clean-versionlib