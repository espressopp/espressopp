=======================
MPI for Python (mpi4py)
=======================

SIAM CSE 2009, Miami
--------------------

* Brian Granger (Cal Poly San Luis Obispo)
* Lisandro Dalcin (UNL, Argentina)

.. class:: center

   http://mpi4py.googlecode.com

   mpi4py@googlegroups.com

.. .. footer:: B. Granger - L. Dalcin, SIAM CSE 2009, Miami


Who are we?
===========

**Lisandro Dalcin**:
    Mastermind, creator and lead developer of MPI For Python.  Also
    developer of petsc4py and slepc4py.  Computational fluid dynamics.

**Brian Granger (speaker)**:
    I help out on a very part time basis with documentation, building,
    testing, design ideas, etc.  Theoretical atomic and molecular
    physics.

Lisandro couldn't attend this meeting and asked me to present his
talk.

What is MPI for Python?
=======================

.. class:: incremental

* Full-featured Python bindings for MPI.

* API based on the standard MPI-2 C++ bindings.

* Almost all MPI calls are supported.

  * targeted to MPI-2 implementations.
    
  * also works with MPI-1 implementations.


Implementation
==============

Implemented with Cython http://www.cython.org

.. class:: incremental

* Code base far easier to write, maintain, and extend.

* Faster than other solutions (mixed Python and C codes).

* A *pythonic* API that runs at C speed !


Portability
===========

.. class:: incremental

* Tested on all major platforms (Linux, Mac OS X, Windows).

* Works with the open-source MPI's (MPICH2, Open MPI, MPICH1, LAM).

* Should work with vendor-provided MPI's (HP, IBM, SGI).

* Works on Python 2.3 to 3.0 (Cython is just great!).


Interoperability
================

Good support for wrapping other MPI-based codes.

.. class:: incremental

* You can use Cython (``cimport`` statement).

* You can use SWIG (*typemaps* provided).

* You can use F2Py (``py2f()``/``f2py()`` methods).

* You can use hand-written C (C-API provided).

mpi4py will allow you to use virtually any MPI based C/C++/Fortran
code from Python.


Features
========

.. class:: incremental
  
* Classical MPI-1 Point-to-Point.

  + blocking (send/recv)

  + non-blocking (isend/irecv, test/wait).

* Classical MPI-1 and Extended MPI-2 Collectives.

* Dynamic Process Management (spawn, accept/connect).

* Parallel I/O (files, read/write).

* One-sided (windows, get/put/accumulate).


Features (cont.)
================

.. class:: incremental

* Communication of general Python objects (pickle).

  .. class:: incremental

  * very convenient, as general as pickle can be.

  * can be slow for large data (CPU and memory consuming).

* Communication of array data (Python's buffer interface).
  
  .. class:: incremental

  * MPI datatypes have to be explicitly specified.

  * very fast, almost C speed (for messages above 5-10 kB).


Features (cont.)
================

* Integration with IPython

  .. sourcecode:: bash

     $ ipcluster mpirun -n 16 --mpi=mpi4py

  .. class:: incremental
  
  * enables MPI applications to be used *interactively*.



Point-to-Point, Python objects
==============================

.. sourcecode:: python

   from mpi4py import MPI

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()

   if rank == 0:
      data = {'a': 7, 'b': 3.14}
      comm.send(data, dest=1, tag=11)
   elif rank == 1:
      data = comm.recv(source=0, tag=11)


Point-to-Point, (NumPy) array data
==================================
.. sourcecode:: python

   from mpi4py import MPI
   import numpy

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()

   if rank == 0:
      data = numpy.arange(1000, dtype='i')
      comm.Send([data, MPI.INT], dest=1, tag=77)
   elif rank == 1:
      data = numpy.empty(1000, dtype='i')
      comm.Recv([data, MPI.INT], source=0, tag=77)


Broadcasting a Python dictionary
================================

.. sourcecode:: python

   from mpi4py import MPI

   comm = MPI.COMM_WORLD
   rank = comm.Get_rank()

   if rank == 0:
      data = {'key1' : [7, 2.72, 2+3j],
              'key2' : ( 'abc', 'xyz')}
   else:
      data = None
   data = comm.bcast(data, root=0)


Parallel matrix-vector product
==============================

.. sourcecode:: python

   from mpi4py import MPI
   import numpy

   def matvec(comm, A, x):
       m = A.shape[0] # local rows
       p = comm.Get_size()
       xg = numpy.zeros(m*p, dtype='d')
       comm.Allgather([x,  MPI.DOUBLE],
                      [xg, MPI.DOUBLE])
       y = numpy.dot(A, xg)
       return y


Compute Pi (master side)
========================

.. class:: tiny
.. sourcecode:: python

   #! /usr/bin/python
   from mpi4py import MPI
   import numpy
   import sys

   comm = MPI.COMM_SELF.Spawn(sys.executable,
                              args=['cpi.py'],
                              maxprocs=5)

   N = numpy.array(100, 'i')
   comm.Bcast([N, MPI.INT], root=MPI.ROOT)
   PI = numpy.array(0.0, 'd')
   comm.Reduce(None, [PI, MPI.DOUBLE],
               op=MPI.SUM, root=MPI.ROOT)
   print(PI)

   comm.Disconnect()


Compute Pi (worker side)
=================================

.. class:: tiny
.. sourcecode:: python

   #! /usr/bin/python
   from mpi4py import MPI
   import numpy

   comm = MPI.Comm.Get_parent()
   size = comm.Get_size()
   rank = comm.Get_rank()

   N = numpy.array(0, dtype='i')
   comm.Bcast([N, MPI.INT], root=0)
   h = 1.0 / N; s = 0.0
   for i in range(rank, N, size):
       x = h * (i + 0.5)
       s += 4.0 / (1.0 + x**2)
   PI = numpy.array(s * h, dtype='d')
   comm.Reduce([PI, MPI.DOUBLE], None,
               op=MPI.SUM, root=0)

   comm.Disconnect()


Throughput, *Sendrecv* exchange
===============================

.. image:: sendrecv1.png
   :align: center 
.. :scale: 60


Overhead, *Sendrecv* exchange
=============================

.. image:: sendrecv2.png
   :align: center 
.. :scale: 60


Wrapping with SWIG
==================

+---------------------------------------+-----------------------------------------+
| .. class:: tiny                       | .. class:: tiny                         |
| .. sourcecode:: none                  | .. sourcecode:: none                    |
|                                       |                                         |
|    // file: helloworld.i              |    /* file: helloworld.c */             |
|    %module helloworld                 |    void sayhello(MPI_Comm comm)         |
|    %{                                 |    {                                    |      
|    #include <mpi.h>                   |      int size, rank;                    |
|    #include "helloworld.c"            |      MPI_Comm_size(comm, &size);        |
|    }%                                 |      MPI_Comm_rank(comm, &rank);        |
|                                       |      printf("Hello, World! "            |
|    %include mpi4py/mpi4py.i           |             "I am process %d of %d.\n", |
|    %mpi4py_typemap(Comm, MPI_Comm);   |             rank, size);                |
|    void sayhello(MPI_Comm comm);      |    }                                    |
+---------------------------------------+-----------------------------------------+
| .. class:: tiny                                                                 |
| .. sourcecode:: python                                                          |
|                                                                                 |
|    >>> from mpi4py import MPI                                                   |
|    >>> import helloworld                                                        |
|    >>> helloworld.sayhello(MPI.COMM_WORLD)                                      |
|    Hello, World! I am process 0 of 1.                                           |
+---------------------------------------------------------------------------------+  


Wrapping with F2Py
==================

+---------------------------------------------------------------------------------+
| .. class:: tiny                                                                 |
| .. sourcecode:: none                                                            |
|                                                                                 |
|    !file: helloworld.f90                                                        |
|    subroutine sayhello(comm)                                                    |
|      use mpi                                                                    |
|      implicit none                                                              |
|      integer :: comm, rank, size, ierr                                          |
|      call MPI_Comm_size(comm, size, ierr)                                       |
|      call MPI_Comm_rank(comm, rank, ierr)                                       |
|      print *, 'Hello, World! I am process ',rank,' of ',size,'.'                |
|    end subroutine sayhello                                                      |
+---------------------------------------------------------------------------------+
| .. class:: tiny                                                                 |
| .. sourcecode:: python                                                          |
|                                                                                 |
|    >>> from mpi4py import MPI                                                   |
|    >>> import helloworld                                                        |
|    >>> fcomm = MPI.COMM_WORLD.py2f()                                            |
|    >>> helloworld.sayhello(fcomm)                                               |
|    Hello, World! I am process 0 of 1.                                           |
+---------------------------------------------------------------------------------+


Disclaimer
==========

Linsandro says, "it is really hard to test all the possibilities..."

* Python versions (2.3 to 3.0)
* MPI implementations (MPI-1/2, open-source, vendor-provided)
* Compilers (GNU, Intel, PathScale), operating systems, batch systems.

Brian says, "true, but..."

* Any parallel, MPI-based code has these problems.
* The mpi4py test suite of mpi4py is *really good*.


Conclusions
===========

* Python is a great language for HPC.
* In addition to mpi4py and IPython there are also: petsc4py,
  slepc4py, pytrilinos.
* Great glue language for mixed language parallel codes.

Do not hesitate to ask for help ...

* Mailing List mpi4py@googlegroups.com

* Lisandro Dalcin dalcinl@gmail.com
