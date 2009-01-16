## THE SCRIPTS
python_SCRIPTS=$(python_scripts)

python_scripts=hello.py 

# python script handling
include $(top_srcdir)/build-aux/macros/axes_python.am

# clean the __init__.py-files from the build dir
MOSTLYCLEANFILES=__init__.py __init__.pyc __init__.install.py

# rename __init__.build.py to __init__.py in the build dir
# copy __init__.install.py into the build dir
all-local: python-all
	sed -e "s,[@]cppdir[@],'$(abs_top_builddir)/src/cpp',g" \
	$(srcdir)/__init__.build.py > __init__.py
	cp $(srcdir)/__init__.install.py .

# install __init__.install.py into pkglibdir, renaming it to __init__.py
install-exec-local:
	$(install_sh_SCRIPT) __init__.install.py $(DESTDIR)$(pythondir)/__init__.py

# remove __init__.py and __init__.pyc from pkglib
uninstall-local:
	-rm $(DESTDIR)$(pythondir)/__init__.py
	-rm $(DESTDIR)$(pythondir)/__init__.pyc



