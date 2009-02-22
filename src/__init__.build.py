# Init file for the ESPResSo package
# This file is used in the build directory only, it will not be
# installed. For the installed version, please have a look at the file
# src/python/__init__.install.py in the ESPResSo sources.
# 
# This file is necessary to make the directory a package, and to
# be able to import the espresso C++ library _espresso.so from the cpp
# source directory.

# Import the shared library from the src dir.
import os, sys
pathname=os.path.join(@srcdir@, '.libs')
sys.path.append(pathname)
