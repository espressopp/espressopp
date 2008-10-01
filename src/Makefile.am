SUBDIRS = cpp

# - make the python script part available in the build tree's src,
#   that is, link it from the source tree into the build tree
# - in the build tree's src, make a link from cpp to espresso,
#   so that importing "espresso" from the build tree's src works
all-local:
	if test ! -e python || test -h python; then \
		rm -f python; \
		$(LN_S) $(srcdir)/python python; \
	else \
		echo "python exists and is not a link, please remove" 1>&2; \
		exit -1; \
	fi
	if test ! -e espresso || test -h espresso; then \
		rm -f espresso; \
		$(LN_S) cpp espresso; \
	else \
		echo "espresso exists and is not a link, please remove" 1>&2; \
		exit -1; \
	fi
clean-local:
	if test -h python; then \
		rm -f python; \
	fi
	if test -h espresso; then \
		rm -f espresso; \
	fi

