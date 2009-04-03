DIST_SUBDIRS=unittest

subpkgincludedir=$(pkgincludedir)/esutil

subpkginclude_HEADERS = \
	SmallVector.hpp \
	DerivableSmallVector.hpp \
	Vector3D.hpp \
	Identifier.hpp \
	TupleVector.hpp \
	Timer.hpp \
	virtual_functional.hpp \
        Logger.hpp

noinst_LTLIBRARIES = libesutil.la \
                     libgenlogger.la \
                     libpylogger.la

libesutil_la_SOURCES = \
	TupleVector.cpp \
	Timer.cpp 

libgenlogger_la_SOURCES = GenLogger.cpp GenLogger.hpp

libpylogger_la_SOURCES = PyLogger.cpp PyLogger.hpp

