SUBDIRS = external modules

EXTRA_DIST = tests.txt

.PHONY: check-extended tests.txt

all-local: tests.txt

# block recursive make check, since make check here will in fact
# recursively check all subdirectories
check-recursive:

# make a link here to tests.txt, so that runtest.sh can be run directly
# without specifying the location of the test files. According to the gnu
# standards, we should not encode the source tree directory into runtests.sh
tests.txt:
	test -e tests.txt || $(LN_S) $(srcdir)/tests.txt .

localtest= $(top_builddir)/bin/runtests.sh
check: all tests.txt
	$(localtest)
check-extended: all tests.txt
	$(localtest) -notag extended
