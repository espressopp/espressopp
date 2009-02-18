subpkgincludedir=$(pkgincludedir)/pmi
subpkginclude_HEADERS = \
        types.hpp \
        transmit.hpp \
        internal_func.hpp \
        functions.hpp \
        exceptions.hpp \
        ParallelMethod.hpp \
        ParallelClass.hpp \
        pmi.hpp

noinst_HEADERS = \
        worker_internal.hpp

# AM_CPPFLAGS = -DLOG4ESPP_LEVEL_OFF

lib_LTLIBRARIES = libpmi.la
libpmi_la_SOURCES = \
        types.cpp \
        internal_func.cpp \
        functions.cpp \
        exceptions.cpp \
        ParallelClass.cpp \
        ParallelMethod.cpp \
        worker_internal.cpp \
        mpitransmit/transmit.cpp
