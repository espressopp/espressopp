SUBDIRS = external modules

EXTRA_DIST = tests.cfg

.PHONY: check-extended

all-local: tests.cfg

# block recursive make check, since make check here will in fact
# recursively check all subdirectories
# TODO: However, it will not BUILD the required files
check-recursive:

check: tests.cfg runtests.sh
	test -e tests.cfg || $(LN_S) $(srcdir)/tests.cfg .
	./runtests.sh
check-extended: tests.cfg runtests.sh
	test -e tests.cfg || $(LN_S) $(srcdir)/tests.cfg .
	./runtests.sh -notag extended
