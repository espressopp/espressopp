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

#################################################################
# In the build tree's src, make a link from src to espresso,
# so that importing "espresso" from the build dir works
all-local:
	if test ! -e espresso || test -h espresso; then \
		rm -f espresso; \
		$(LN_S) src espresso; \
	else \
		echo "espresso exists and is not a link, please remove" 1>&2; \
		exit -1; \
	fi
clean-local:
	if test -h espresso; then \
		rm -f espresso; \
	fi

