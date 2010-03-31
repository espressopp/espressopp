================================
MPI for Python (mpi4py) *part I*
================================

.. raw:: html   

   <p class="center">Eilif Muller (EPFL, Switzerland) and<br>Lisandro Dalcin (UNL, Argentina)</p>


.. class:: center

   http://mpi4py.googlecode.com

.. footer:: 
   .. raw:: html

      <h1>MPI for Python</h1>
      <h2>Muller & Dalcin, CodeJam3, Freiburg</h2>


.. header:: 
   
   .. raw:: html   
   
	<img class="scale" align="right" width="140" height="95" src="logo.png" alt="mpi4py" title="mpi4py" />


Who are we?
===========

**Lisandro Dalcin**:
    Mastermind, creator and lead developer of MPI For Python.  Also
    developer of petsc4py and slepc4py.  Computational fluid dynamics.

**Eilif Muller (speaker)**:
    After evaluating other Python MPI packages (pypar, pyMPI, etc.), I use mpi4py with IPython 
    since 2004 for developing various parallelized programs in Python.

.. class:: nb

   Lisandro couldn't accept our invitation to attend and asked me to present in his place.

Multi-core crisis
=================

.. raw:: html   
   
    <div align="center" class="align-center"><img class="scale" width="700" height="520" src="clockspeeds.jpg" alt="clockspeeds" title="clockspeeds" /></div>


The free lunch is over
======================

.. class:: center blue huge  

   “Concurrency is the next major revolution in how we write software [after OOP].”

.. list-table::
  :class: borderless

  * - Herb Sutter, *The Free Lunch is Over: 
      A Fundamental Turn Towards Concurrency in Software*, 
      Dr.Dobb's Journal, 30(3) March 2005.

Three types of Concurrency
==========================

.. list-table::

  * - **SMP**
    - **Message passing**
    - **GPU/Stream**

  * - Shared mem., Threads

    - MPI standard, Sockets

    - GPU, Cell

  * - <8 threads, or $

    - Linux clusters: 1000's of processes over network

    - Stream kernels over arrays, upto 320 hardware threads @ float5 SIMD (ATI 5870)

  * - threads, multiprocessing

    - mpi4py, sockets

    - PyCUDA, PyOpenCL



Message passing Concurrency
===========================


* Multi-process execution facilities

  .. sourcecode:: bash

     $ mpiexec -n 16 python helloworld.py


* API for inter-process message exchanges
  
  * eg. Basic P2P: Send (emitter), Recv (consumer)

What is MPI for Python?
=======================

.. class:: incremental

* A wrapper for widely used MPI (MPICH2, OpenMPI, LAM/MPI)

  * MPI supported by wide range of vendors, hardware, languages

* API based on the standard MPI-2 C++ bindings.

* Almost all MPI calls are supported.

  * targeted to MPI-2 implementations.
    
  * also works with MPI-1 implementations.


Basic stuff
===========

.. sourcecode:: python

   from mpi4py import MPI

* Communicator = Comm
  
  * Manages processes and communication between them

* MPI.COMM_WORLD

  * all processes defined at exec.time

* Comm.size, Comm.rank

Basic stuff (cont.)
===================

.. sourcecode:: python

   from mpi4py import MPI
   comm = MPI.COMM_WORLD
   print "Hello from %s, %d of %d"\
     % (MPI.Get_processor_name(),
        comm.rank, comm.size)

→ test.py

:: 

   $ mpiexec -n 2 python test.py
   Hello from rucola, 0 of 2
   Hello from rucola, 1 of 2


Point-to-Point: Python objects
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


P2P: (NumPy) array data
=======================

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


Throughput, *Sendrecv* exchange
===============================

.. raw:: html

   <div align="center" class="align-center"><img class="scale" width="662" height="530" src="sendrecv1.png" alt="sendrecv1.png" title="sendrecv1.png" /></div>   



Overhead, *Sendrecv* exchange
=============================

.. raw:: html

   <div align="center" class="align-center"><img class="scale" width="662" height="530" src="sendrecv2.png" alt="sendrecv2.png" title="sendrecv2.png" /></div>   



Non-blocking P2P
================

.. sourcecode:: python

  if rank == 0:
     data = numpy.arange(1000, dtype='i')
     req = comm.Isend([data, MPI.INT],dest=1,...)
  elif rank == 1:
     data = numpy.empty(1000, dtype='i')
     req1 = comm.Irecv([data, MPI.INT],source=0,...)

.. class:: blue

  < do something, even another Irecv, etc. >

.. sourcecode:: python

  if rank == 1:
     status = [MPI.Status(), ... ]
     MPI.Request.Waitall([req1, ...], status)

Persitant P2P
=============

Store messaging paramters as a Prequest to be used in a loop:

.. sourcecode:: python

   request = comm.Recv_init([msg,MPI.INT],
			    partner_rank)
   for i in xrange(10):

       MPI.Prequest.Startall(request)

       do_something()

       MPI.Request.Waitall([request])

       do_something_with(msg)


Collective Messages
===================
 
Involve the whole Comm

* ::

   Scatter
   - Spread a sequence over processes
   
* ::

   Gather
   - Collect a sequence scattered over 
     processes

* ::

   Broadcast
   - Send a message to all processes

* ::
   
   Barrier 
   - block till all processes arrive

Scatter & Gather
================

.. sourcecode:: python

   N = 100
   assert(N%com.size==0)
   if com.rank==0: msg = numpy.arange(N,dtype=float)
   else: msg = None

   dest = numpy.empty(N/com.size, dtype=float)
   ans = numpy.empty(com.size, dtype=float)

   com.Scatter([msg,MPI.DOUBLE],
               [dest,MPI.DOUBLE],root=0)
   mysum = numpy.sum(dest)

   com.Gather([mysum, MPI.DOUBLE], 
   	      [ans,MPI.DOUBLE], root=0)
   if com.rank==0: print ans

ans → [1225. 3725.]


Array data buffer notation
==========================

.. class:: blue

  Basic: [buf, MPI datatype]

.. sourcecode:: python

   a = numpy.empty(10,dtype=float)
   comm.Send([a, MPI.DOUBLE], dest=1, tag=77)    

.. class:: blue

  Vector collectives: [buf, count, displ, MPI datatype]

.. sourcecode:: python

   comm.Scatterv([msg,counts,None,MPI.DOUBLE],
		 [b,MPI.DOUBLE])
   comm.Allgatherv([b,MPI.DOUBLE],
                   [c,counts,None,MPI.DOUBLE])

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

* You can use boost.

* You can use SWIG (*typemaps* provided).

* You can use F2Py (``py2f()``/``f2py()`` methods).

* You can use hand-written C (C-API provided).

mpi4py will allow you to use virtually any MPI based C/C++/Fortran
code from Python.

.. class:: blue

  More on this in part II

Features Summary
================

.. class:: incremental
  
* Classical MPI-1 Point-to-Point.

  + blocking (send/recv)

  + non-blocking (isend/irecv, test/wait).

* Classical MPI-1 and Extended MPI-2 Collectives.

.. class:: blue

  More on these in part II:

* Dynamic Process Management (spawn, accept/connect).

* Parallel I/O (files, read/write).

* One-sided (windows, get/put/accumulate).


Features Summary (cont.)
========================

.. class:: incremental

* Communication of general Python objects (pickle).

  .. class:: incremental

  * very convenient, as general as pickle can be.

  * can be slow for large data (CPU and memory consuming).

* Communication of array data (Python's buffer interface).
  
  .. class:: incremental

  * MPI datatypes have to be explicitly specified.

  * very fast, almost C speed (for messages above 5-10 kB).


Features Summary (cont.)
========================

* Integration with IPython

  .. sourcecode:: bash

     $ ipcluster mpirun -n 16 --mpi=mpi4py

  .. class:: incremental
  
  * enables MPI applications to be used *interactively*.



Disclaimer
==========

Linsandro says, "it is really hard to test all the possibilities..."

* Python versions (2.3 to 3.0)
* MPI implementations (MPI-1/2, open-source, vendor-provided)
* Compilers (GNU, Intel, PathScale), operating systems, batch systems.

But ...

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

.. class:: blue

  More on DPM, Wrapping, Reduce Ops, Parallel I/O in part II

