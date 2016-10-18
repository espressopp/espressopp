#!/usr/bin/env python
import os
os.environ['MPE_LOGFILE_PREFIX'] = 'ring'
import mpi4py
mpi4py.profile('mpe')

from mpi4py import MPI
from array import array

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

src  = rank-1
dest = rank+1
if rank == 0:
    src = size-1
if rank == size-1:
    dest = 0

try:
    from numpy import zeros
    a1 = zeros(1000000, 'd')
    a2 = zeros(1000000, 'd')
except ImportError:
    from array import array
    a1 = array('d', [0]*1000); a1 *= 1000
    a2 = array('d', [0]*1000); a2 *= 1000

comm.Sendrecv(sendbuf=a1, recvbuf=a2,
              source=src, dest=dest)

MPI.Request.Waitall([
    comm.Isend(a1, dest=dest),
    comm.Irecv(a2, source=src),
    ])
