# Init file for the ESPResSo package
# This file is used in the build directory only, it wil not be
# installed. For the installed version, please have a look at the file
# src/python/__init__.install.py in the ESPResSo sources.
# 
# This file is necessary to make the directory a package, and to
# import the espresso C++ library _escpp.so from the cpp source
# directory.

# Import the libtool library from the cpp dir.
# This code was adapted from ltihooks.py by James Henstridge in the
# pygtk package.
import os, imp

# get the name of the actual dynlib from the .la file
lafile=os.path.join(@cppdir@,'_escpp.la')
fp = open(lafile, 'r');
dlname = ''
line = fp.readline()
while line:
    if len(line) > 7 and line[:7] == 'dlname=':
        dlname = line[8:-2]
    line = fp.readline()
fp.close()
if dlname:
    filename = os.path.join(@cppdir@,
                            '.libs', dlname)

# now load the dynlib
imp.load_dynamic("_escpp", filename)
