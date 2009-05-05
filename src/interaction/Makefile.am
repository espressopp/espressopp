subpkgincludedir=$(pkgincludedir)/interaction
subpkginclude_HEADERS = \
	bindings.hpp \
	Interaction.hpp \
	LennardJones.hpp \
	FENE.hpp

noinst_LTLIBRARIES = libinteraction.la
libinteraction_la_SOURCES = \
	bindings.cpp \
	Interaction.cpp \
	LennardJones.cpp \
        FENE.cpp

include $(top_srcdir)/build-aux/macros/axes_python.am
python_SCRIPTS = __init__.py
