Installation of ESPResSo++
==========================

The first step in the installation of ESPResSo++ is to download it from the
following location:

     https://www.espressopp-pp.de/Download/espressopp_1_3_1.tgz

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
Alternatively, if cmake . complains on missing BOOST or MPI4PY libraries and you had not
installed them, you can try

# cmake . -DEXTERNAL_BOOST=OFF -DEXTERNAL_MPI4PY=OFF

In this case, ESPResSo++ will try to use internal Boost and mpi4py libraries.

After successfully building all the Makefiles you should build ESPResSo++ with:

# make
(This will take several minutes)

Before beeing able to use the espressopp  module in Python you need to source the ESPRC file:

# source ESPRC
(This sets all corresponding environment variables to point to the module, e.g. PYTHONPATH)
You have to source this file every time you want to work with espressopp. It would advisable to
e.g. source the file in your .bashrc file ( "source <path_to_espressopp>/ESPRC" )

In order to use matplotlib.pyplot for graphical output get the open source code from:

  http://sourceforge.net/projects/matplotlib

and follow the installation instructions of your distribution.

