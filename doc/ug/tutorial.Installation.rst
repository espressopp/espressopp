Installation of ESPResSo++
==========================

The first step in the installation of ESPResSo++ is to download it from the
following location:

     https://www.espresso-pp.de/Download/espressopp_1_3_1.tgz

On the command line type:

# tar -xzf espressopp_1_3_1.tgz

This will create a subdirectory espressopp-1.3.1

Enter this subdirectory

# cd espressopp-1.3.1

Create the Makefiles using the cmake command. If you don't have it yet, you have to
install it first. It is available for all major Linux distributions and also for Mac OS X.
(ubuntu,debian: "apt-get install cmake" or get it from http://www.cmake.org )

# cmake .
(the space and dot after *cmake* are necessary)

If cmake doesn't finish sucessfully (e.g. it didn't find all the libraries) you can
tell cmake manually, where to find them by typing:

# ccmake .

This will open an interactive page where all configuration information can be specified.
After successfully building all the Makefiles you should build ESPResSo++ with:

# make
(This will take several minutes)

In order to use matplotlib.pyplot for graphical output get the open source code from:

  http://sourceforge.net/projects/matplotlib

and follow the installation instructions of your distribution.

