SUBDIRS=unittest

subpkgincludedir=$(pkgincludedir)/esutil

subpkginclude_HEADERS = \
	bindings.hpp \
	Vector3D.hpp \
	SmallVector.hpp \
	DerivableSmallVector.hpp \
	Identifier.hpp \
	TupleVector.hpp \
	Timer.hpp \
	virtual_functional.hpp

noinst_LTLIBRARIES = libesutil.la

libesutil_la_SOURCES = \
        bindings.cpp \
	TupleVector.cpp \
	Timer.cpp 
