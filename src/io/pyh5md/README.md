pyh5md : Read and write H5MD files
==================================

Copyright © 2012-2016 Pierre de Buyl and contributors  
*License:* BSD

pyh5md is a library to read and write easily H5MD files. [H5MD][] is a file
format specification based on [HDF5][] to store molecular data. pyh5md is built
on top of [h5py][], the HDF5 for Python library by Andrew Colette.

[H5MD]: http://nongnu.org/h5md/
[HDF5]: http://www.hdfgroup.org/HDF5/
[h5py]: http://h5py.org/

Install
-------

    python setup.py install --user

or

    pip install --user .

installs pyh5md for the current user

Examples
--------

Once pyh5md is installed:

    cd examples
    python random_walk_1d.py

executes an example program that generates the H5MD file `walk_1d.h5`.
