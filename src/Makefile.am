SUBDIRS = cpp python

# in the build tree's src, make a link from python to espresso,
# so that importing "espresso" from the build tree's src works
all-local:
	if test ! -e espresso || test -h espresso; then \
		rm -f espresso; \
		$(LN_S) python espresso; \
	else \
		echo "espresso exists and is not a link, please remove" 1>&2; \
		exit -1; \
	fi
clean-local:
	if test -h espresso; then \
		rm -f espresso; \
	fi
