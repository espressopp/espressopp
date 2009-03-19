ACLOCAL_AMFLAGS = -I build-aux/macros
SUBDIRS = contrib src test doc
EXTRA_DIST = autogen.sh 

#################################################################
# Documentation
#################################################################
.PHONY: FORCE doc ug dg
doc ug dg: FORCE
	cd doc; $(MAKE) --print-directory $@
FORCE:

distclean-local:
	-rm -f config.summary
# remove the espresso and examples links created in configure
	if test -h espresso; then \
		rm -f espresso; \
	fi

install-exec-local:
# install espresso setup scripts
	$(install_sh_SCRIPT) espresso_setup.sh.install $(DESTDIR)$(pkglibdir)/espresso_setup.sh
	$(install_sh_SCRIPT) espresso_setup.py.install $(DESTDIR)$(pkglibdir)/espresso_setup.py

uninstall-local:
	rm -f $(DESTDIR)$(pkglibdir)/espresso_setup.sh
	rm -f $(DESTDIR)$(pkglibdir)/espresso_setup.py

