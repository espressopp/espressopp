libtool: $(LIBTOOL_DEPS)
	$(SHELL) ./config.status --recheck

lib_LTLIBRARIES = libboost_mpi.la
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
