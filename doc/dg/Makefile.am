EXTRA_DIST=doxygen.cfg

.PHONY: doc dg
doc: dg
dg:
if HAVE_DOXYGEN
	SRCDIR=$(top_srcdir) PACKAGE="$(PACKAGE_STRING)"\
	  $(DOXYGEN) $(srcdir)/doxygen.cfg \
	  | grep -ve "^\(Generating\|Parsing\|Preprocessing\)"
	@echo "***************************************************************************"
	@echo "* The Developer's Guide is now accessible at"
	@echo "*  file://`pwd`/html/index.html"
	@echo "***************************************************************************"
else
	@echo "doxygen was not found in your PATH."
	@echo "Can't build developer's guide without doxygen."
	@echo "Install doxygen and rerun configure."

endif

clean-local:
	-rm -rf html
