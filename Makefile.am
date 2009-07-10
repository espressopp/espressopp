ACLOCAL_AMFLAGS = -I build-aux/macros
SUBDIRS = contrib src test doc
EXTRA_DIST = autogen.sh 

noinst_SCRIPTS = \
	bin/utexec \
	bin/espp_mpiexec.sh \
	espresso_setup.sh \
	espresso_setup.sh.install \
	espresso_setup.py \
	espresso_setup.py.install \
	espresso_setup.csh

#################################################################
# Documentation
#################################################################
.PHONY: FORCE doc ug dg
doc ug dg: FORCE
	cd doc; $(MAKE) $(AM_MAKEFLAGS) --print-directory $@
FORCE:

distclean-local:
	-rm -f config.summary
# remove the espresso and examples links created in configure
	if test -h espresso; then \
		rm -f espresso; \
	fi

install-exec-local: espresso_setup.sh.install espresso_setup.py.install
# install espresso setup scripts
	$(install_sh_SCRIPT) espresso_setup.sh.install $(DESTDIR)$(pkglibdir)/espresso_setup.sh
	$(install_sh_SCRIPT) espresso_setup.py.install $(DESTDIR)$(pkglibdir)/espresso_setup.py

uninstall-local:
	rm -f $(DESTDIR)$(pkglibdir)/espresso_setup.sh
	rm -f $(DESTDIR)$(pkglibdir)/espresso_setup.py

