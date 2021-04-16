Parallel H5MD Output
====================

H5MD_ is a format specification on top of the HDF5_ file format.

.. _H5MD: https://nongnu.org/h5md/
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/

Dump Simulation
---------------

Dumping a simulation with ``espressopp.io.DumpH5MDParallel``

::

    dump_h5md_parallel = espressopp.io.DumpH5MDParallel(system, 'dump.h5')
    dump_h5md_parallel.dump()


Restore Simulation
------------------

Restoring a simulation with ``espressopp.io.RestoreH5MDParallel``

::

    system.storage.removeAllParticles()
    restore_h5md_parallel = espressopp.io.RestoreH5MDParallel(system, 'dump.h5')
    restore_h5md_parallel.restore()
    system.storage.decompose()

**ATTENTION**
 *  No checks for duplicates are performed nor is the particle storage cleared
    before inserting new particles. You might want to remove all particles from
    the simulation before calling ``restore``.
 *  Particles are inserted equally among the processes without obeying subdomains.
    You have to call ``decompose`` to ensure particles are located on the correct
    process after restoration.