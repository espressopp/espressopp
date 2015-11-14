Logging mechanism
=================

.. |espp| replace:: ESPResSo++

|espp| uses Loggers

Logging can be switched on in your python script with the following command:

>>> logging.getLogger("*name of the logger*").setLevel(logging.*Level*)

*Level* is one of the following:

======= =====================================================================
 ERROR   for errors that might still allow the application to continue
 WARN    for potentially harmful situations
 INFO    informational messages highlighting progress
 DEBUG   designates fine-grained informational events
======= =====================================================================

Example:

>>> import espressopp
>>> import logging
>>> logging.getLogger("Storage").setLevel(logging.ERROR)

To log everything (WARNING: this will produce **lots** of output):

>>> logging.getLogger("").setLevel(logging.DEBUG)

The following loggers are currently available:

- Configurations
- Observable
- Velocities
- BC
- Logger
- FixedListComm
- FixedPairList
- FixedQuadrupleList
- FixedTripleList
- FixedTupleList
- Langevin
- MDIntegrator
- AngularPotential
- DihedralPotential
- Interaction
- InterpolationAkima
- InterpolationCubic
- InterpolationLinear
- InterpolationTable
- Potential
- CellListAllPairsIterator
- DomainDecomposition.CellGrid
- DomainDecomposition
- DomainDecomposition.NodeGrid
- Storage
- DomainDecompositionAdress
- StorageAdress
- VerletList
- VerletList

