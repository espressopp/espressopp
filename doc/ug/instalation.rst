.. |espp| replace:: ESPResSo++

Installation
==========================

The first step in the installation of |espp| is to download the latest release from the
following location:

     https://github.com/espressopp/espressopp/releases

On the command line type:

.. parsed-literal::

   tar -xzf espressopp-|version|.tgz

This will create a subdirectory espressopp-|version|

Enter this subdirectory

.. parsed-literal::

   cd espressopp-|version|

Create the Makefiles using the cmake command. If you don't have it yet, you have to
install it first. It is available for all major Linux distributions and also for Mac OS X.
(ubuntu,debian: "apt-get install cmake" or get it from http://www.cmake.org )

.. code-block:: bash

   cmake .

(the space and dot after *cmake* are necessary)

If cmake doesn't finish successfully (e.g. it didn't find all the libraries) you can
tell cmake manually, where to find them by typing:

.. code-block:: bash

   ccmake .

This will open an interactive page where all configuration information can be specified.
Alternatively, if cmake . complains on missing BOOST or MPI4PY libraries and you had not
installed them, you can try

.. code-block:: bash

   cmake . -DEXTERNAL_BOOST=OFF -DEXTERNAL_MPI4PY=OFF

In this case, |espp| will try to use internal Boost and mpi4py libraries.

After successfully building all the Makefiles you should build |espp| with:

.. code-block:: bash

   make

(This will take several minutes)

Before beeing able to use the espressopp  module in Python you need to source the ESPRC file:

.. code-block:: bash

   source ESPRC

(This sets all corresponding environment variables to point to the module, e.g. PYTHONPATH)
You have to source this file every time you want to work with espressopp. It would advisable to
e.g. source the file in your .bashrc file ( "source <path_to_espressopp>/ESPRC" )

In order to use matplotlib.pyplot for graphical output get the open source code from:

  http://sourceforge.net/projects/matplotlib

and follow the installation instructions of your distribution.

