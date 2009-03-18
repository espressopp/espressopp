SUBDIRS=unittest

subpkgincludedir=$(pkgincludedir)/esutil

subpkginclude_HEADERS = \
	bindings.hpp \
	estypes.hpp \
	SmallVector.hpp \
	DerivableSmallVector.hpp \
	Vector3D.hpp \
	Real3D.hpp \
	Identifier.hpp \
	TupleVector.hpp \
	Timer.hpp \
	virtual_functional.hpp

noinst_LTLIBRARIES = libesutil.la

libesutil_la_SOURCES = \
        bindings.cpp \
	Real3D.cpp \
	TupleVector.cpp \
	Timer.cpp 
