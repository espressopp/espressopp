ACLOCAL_AMFLAGS = -I build-aux/macros
SUBDIRS = src test doc
EXTRA_DIST = autogen.sh 

#################################################################
# Documentation
#################################################################
.PHONY: FORCE doc ug dg
doc ug dg: FORCE
	cd doc; $(MAKE) --print-directory $@
FORCE:

