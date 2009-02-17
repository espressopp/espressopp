# Init file for the ESPResSo package
# This file is used in the build directory only, it wil not be
# installed. For the installed version, please have a look at the file
# src/python/__init__.install.py in the ESPResSo sources.
# 
# This file is necessary to make the directory a package, and to
# import the espresso C++ library _escpp.so from the cpp source
# directory.

# Import the shared library from the src dir.
import os
pathname=os.path.join(@srcdir@, '.libs')
__path__.append(pathname)
