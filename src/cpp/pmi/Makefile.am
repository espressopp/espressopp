subpkgincludedir=$(pkgincludedir)/pmi
nobase_subpkginclude_HEADERS = \
	types.hpp \
	transmit.hpp \
	basic_func.hpp \
	exceptions.hpp \
	ParallelMethod.hpp \
	ParallelClass.hpp \
	ParallelObject.hpp \
	controller_func.hpp \
	worker_func.hpp \
	controller.hpp \
	worker.hpp \
	pmi.hpp

noinst_HEADERS = \
	worker_internal.hpp 

noinst_LTLIBRARIES = libpmi.la 
libpmi_la_SOURCES = \
	types.cpp \
	basic_func.cpp \
	exceptions.cpp \
	ParallelClass.cpp \
	ParallelMethod.cpp \
	ParallelObject.cpp \
	controller_func.cpp \
	worker_func.cpp \
	worker_internal.cpp \
	mpitransmit/transmit.cpp


