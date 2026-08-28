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

ESPResSo++ requires Python 3.7 or newer. All required Python packages are listed in `requirements.txt`. You can install them via:

.. code-block:: bash

   pip3 install -r requirements.txt

Create the Makefiles using the cmake command. If you don't have it yet, you have to
install it first. It is available for all major Linux distributions and also for Mac OS X.
(ubuntu,debian: "apt-get install cmake" or get it from http://www.cmake.org )

.. code-block:: bash

   cmake -B builddir .

(the space and dot after *cmake* are necessary)

If cmake doesn't finish successfully (e.g. it didn't find all the libraries) you can
tell cmake manually, where to find them by typing:

.. code-block:: bash

   ccmake -B builddir .

This will open an interactive page where all configuration information can be specified.

After successfully building all the Makefiles you should build |espp| with:

.. code-block:: bash

   cmake --build builddir

(This will take several minutes)

After successfully building |espp| add the build directory to your PYTHONPATH:

.. code-block:: bash

   export PYTHONPATH=$PWD/builddir${PYTHONPATH:+:$PYTHONPATH}

In order to use matplotlib.pyplot for graphical output get the open source code from:

  https://matplotlib.org/

and follow the installation instructions of your distribution.

