SUBDIRS=cpp

# make the python script part available in the build tree
all-local:
	if ! test -e python; then \
		$(LN_S) $(top_srcdir)/execution/python python; \
	fi
clean-local:
	if test -h python; then \
		rm -f python; \
	fi
