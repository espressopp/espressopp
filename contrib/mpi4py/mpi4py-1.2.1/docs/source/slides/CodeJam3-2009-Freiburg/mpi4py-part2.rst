=================================
MPI for Python (mpi4py) *part II*
=================================

.. raw:: html   

   <p class="center">Eilif Muller (EPFL, Switzerland) and<br>Lisandro Dalcin (UNL, Argentina)</p>


.. class:: center

   http://mpi4py.googlecode.com

.. footer:: 
   .. raw:: html

      <h1>MPI for Python II</h1>
      <h2>Muller & Dalcin, CodeJam3, Freiburg</h2>


.. header:: 
   
   .. raw:: html   
   
	<img class="scale" align="right" width="140" height="95" src="logo.png" alt="mpi4py" title="mpi4py" />


Matrix multiplication
=====================

.. sourcecode:: python

  # get a dimension which is multiple of # of cpus
  x = int(numpy.ceil(4000.0/num_cpus)*num_cpus)
  y = 2000
  step = x/num_cpus

  # scatter destination
  domain = numpy.zeros((step,y),dtype='d')
  # gather destination c = dot(a,b)
  c = numpy.zeros((x,x),dtype='d')

  if rank==0:
     # two matrices to multiply
     a = random.uniform(size=(x,y)).astype('d')
     b = random.uniform(size=(y,x)).astype('d')
  else:
     # target buffers
     a = numpy.zeros((x,y),dtype='d')
     b = numpy.zeros((y,x),dtype='d')


Matrix multiplication
=====================

.. sourcecode:: python

  # broadcast b matrix
  COMM.Bcast([b,MPI.DOUBLE])

  COMM.Scatter([a,MPI.DOUBLE],[domain,MPI.DOUBLE])
  ans = numpy.dot(domain,b)    
  COMM.Gather([ans,MPI.DOUBLE],[c,MPI.DOUBLE])


Reduce
======

Perform PROD, SUM, MIN, MAX over comm

.. class:: small

.. sourcecode:: python

   x = numpy.random.uniform(size=2)
   y = numpy.empty(1,dtype=float)
   comm.Reduce([x, MPI.DOUBLE], [y, MPI.DOUBLE]
               op=MPI.MAX, root=0)
   
   print "x=",x, "on %d, y=%f"% (comm.rank, y[0])

::

 x=[ 0.77836931  0.49683394] on 0, y=0.92012392
 x=[ 0.92012392  0.83533206] on 1, y=0.000000
 x=[ 0.28767185  0.29955473] on 2, y=0.000000
 x=[ 0.89894418  0.55617815] on 3, y=0.000000

.. sourcecode:: python

   comm.Allreduce(MPI.IN_PLACE, [x, MPI.DOUBLE],
                  op=MPI.MIN)

::

 x = [ 0.28767185  0.29955473] on 0,1,2,3


Scan
====

.. sourcecode:: python

   x = numpy.array((rank+1),dtype=float)
   y = numpy.empty(1,dtype=float)
   comm.Scan([x, MPI.DOUBLE], [y, MPI.DOUBLE],
             op=MPI.PROD)

   print y, "on %d"% (comm.rank,)

::

 [ 1.] on 0  [ 2. ] on 1
 [ 6.] on 2  [ 24.] on 3

.. sourcecode:: python

   comm.Exscan([x, MPI.DOUBLE], [y, MPI.DOUBLE],
               op=MPI.PROD)

::

 [ 0.] on 0  [ 1.] on 1  
 [ 2.] on 2  [ 6.] on 3




reduce
======

.. sourcecode:: python

   x = [rank]
   y = comm.reduce(x, None, op=MPI.SUM, root=0)

::

 None on 1,2,3
 [0, 1, 2, 3] on 0

scan
====

.. sourcecode:: python

   x = [rank]
   y = comm.scan(x, None, op=MPI.SUM)

::
 
 [0] on 0 
 [0, 1] on 1
 [0, 1, 2] on 2
 [0, 1, 2, 3] on 3

.. sourcecode:: python

   x = [rank]
   y = comm.exscan(x, None, op=MPI.SUM)

::

 None on 0
 [0] on 1
 [0, 1] on 2
 [0, 1, 2] on 3

Custom Reduce Ops
=================

.. sourcecode:: python

    def mysum_py(a, b):
        for i in range(len(a)):
            b[i] = a[i] + b[i]
        return b

    def mysum(ba, bb, dt):
        if dt is None:
            return mysum_py(ba, bb)
        a = numpy.frombuffer(bytes(ba),dtype='i')
        b = numpy.frombuffer(bytes(bb),dtype='i')
        bb[:] = (b+a).data

    myop = MPI.Op.Create(mysum, commute=True)

.. class:: blue

  <stuff>

.. sourcecode:: python

    myop.Free()

Custom Reduce Ops (cont)
========================

.. sourcecode:: python

    arr = numpy.array([rank],dtype='i')
    dest = numpy.empty(1,dtype='i')

    comm.Reduce([arr,MPI.INT], [dest,MPI.INT], 
    	        op=myop, root=0)
    y = comm.reduce([rank], None, op=myop)

    myop.Free()

::

 On 0: y = [6], dest = array([6])


Dynamic Process Management
==========================

* Process can spawn a new Comm

* New feature of MPI-2 standard

* Possible outside of an mpiexec execution model

.. sourcecode:: python

   MPI.COMM_SELF.Spawn(sys.executable,
                       args=['worker.py'],
                       maxprocs=5)


Compute Pi (master side)
========================

.. class:: small
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

.. class:: small
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


Parallel I/O
============

.. sourcecode:: python

   comm = MPI.COMM_WORLD
   atype = MPI.DOUBLE
   mode = MPI.MODE_WRONLY | MPI.MODE_CREATE

   m, M = 1,1*comm.size
   A = rank*numpy.ones(m*M, dtype='d').reshape(m,M)

   sizes = [M, M] # global array shape
   subsizes = [m, M] # local subarray shape
   starts = [m*comm.rank, 0] # start of section here
   mktype = atype.Create_subarray # constructor
   view = mktype(sizes,subsizes,starts,MPI.ORDER_C)
   fh = MPI.File.Open(comm, 'datafile', mode )
   fh.Set_view(etype=atype, filetype=view)

   fh.Write_all([A, atype])
   fh.Close()
   view.Free()

Created 'datafile'

Parallel I/O (cont)
===================

.. sourcecode:: python

   f = file('datafile','rb')
   a = numpy.frombuffer(f.read(),dtype=float)

a = [ 0.  0.  0.  1.  1.  1.  2.  2.  2.]

* pytables (hdf5) parallel I/O support?



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
|                                                                                 |
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
|                                                                                 |
| .. sourcecode:: python                                                          |
|                                                                                 |
|    >>> from mpi4py import MPI                                                   |
|    >>> import helloworld                                                        |
|    >>> fcomm = MPI.COMM_WORLD.py2f()                                            |
|    >>> helloworld.sayhello(fcomm)                                               |
|    Hello, World! I am process 0 of 1.                                           |
+---------------------------------------------------------------------------------+

Wrapping with boost
===================

.. sourcecode:: c++

   static void hw_sayhello(object py_comm)
   {
      PyObject* py_obj = py_comm.ptr();
      MPI_Comm *comm_p = PyMPIComm_Get(py_obj);
      if (comm_p == NULL) throw_error_already_set();
      sayhello(comm_p);
   }

   BOOST_PYTHON_MODULE(helloworld)
   {
      if (import_mpi4py() < 0) return; /* Python 2.X */

      def("sayhello", hw_sayhello);
   }

.. sourcecode:: python

   import helloworld
   helloworld.sayhello(comm)


Conclusions
===========

* Python is a great language for HPC.
* In addition to mpi4py and IPython there are also: petsc4py,
  slepc4py, pytrilinos.
* Great glue language for mixed language parallel codes.

Do not hesitate to ask for help ...

* Mailing List mpi4py@googlegroups.com

* Lisandro Dalcin dalcinl@gmail.com
