#  Copyright (c) 2015
#      Jakub Krajniak (jkrajniak at gmail.com)
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


r"""
*********************************************
**DumpPairs** - IO Object
*********************************************

* `dump()`

  Properties

* `h5md_file`
  HDF5 file object.
"""

from espressopp.esutil import cxxinit
from espressopp import pmi
from mpi4py import MPI

from espressopp.ParticleAccess import *
from _espressopp import io_DumpTopology

import collections
import pyh5md
import numpy as np


class DumpTopologyLocal(ParticleAccessLocal, io_DumpTopology):
    def __init__(self, system, integrator, h5md_file):
        if pmi.workerIsActive():
            cxxinit(self, io_DumpTopology, system, integrator)
            self.h5md_file = h5md_file
            self.tuple_index = 0
            self.tuple_data = {}
            if 'connectivity' not in self.h5md_file.file.f:
                self.h5md_file.file.f.create_group('connectivity')
            self.chunk_size = 128
            self.dt = integrator.dt

    def dump(self):
        if pmi.workerIsActive():
            self.cxxclass.dump(self)

    def observe_tuple(self, fpl, name, particle_group='atoms'):
        if pmi.workerIsActive():
            self.cxxclass.observe_tuple(self, fpl)
            g = pyh5md.TimeData(
                self.h5md_file.file.f['/connectivity'],
                name,
                shape=(self.chunk_size, 2),
                dtype=np.int,
                fillvalue=-1)
            g.attrs['particle_group'] = particle_group
            self.tuple_data[self.tuple_index] = g
            self.tuple_index += 1

    def update(self):
        if pmi.workerIsActive():
            raw_data = self.cxxclass.get_data(self)
            step_data = collections.defaultdict(dict)
            max_size = 0
            while raw_data != []:
                step = raw_data.pop()
                fpl_idx = raw_data.pop()
                fpl_size = raw_data.pop()
                data = []
                for _ in range(fpl_size):
                    b1 = raw_data.pop()
                    b2 = raw_data.pop()
                    data.append((b1, b2))
                max_size = max(len(data), max_size)
                step_data[step][fpl_idx] = np.array(data, dtype=np.int)
            MPI.COMM_WORLD.Barrier()
            NMaxLocal = np.array(max_size, 'i')
            NMaxGlobal = np.array(0, 'i')
            MPI.COMM_WORLD.Allreduce(NMaxLocal, NMaxGlobal, op=MPI.MAX)
            cpu_size = ((NMaxGlobal // self.chunk_size)+1)*self.chunk_size
            total_size = MPI.COMM_WORLD.size*cpu_size
            idx_0 = MPI.COMM_WORLD.rank*cpu_size
            idx_1 = idx_0 + NMaxLocal
            for step in sorted(step_data):
                for fpl_idx in step_data[step]:
                    data = step_data[step][fpl_idx]
                    if data.shape[0] < NMaxLocal:
                        old_shape = data.shape
                        if len(data):
                            linear_data = data.reshape((old_shape[0]*old_shape[1],))
                        else:
                            linear_data = data
                        data = np.pad(linear_data, (0, 2*NMaxLocal-linear_data.shape[0]),
                                      'constant', constant_values=-1)
                        data = data.reshape((NMaxLocal, 2))
                    ds = self.tuple_data[fpl_idx]
                    if total_size > ds.value.shape[1]:
                        ds.value.resize(total_size, axis=1)
                    ds.append(data, step, step*self.dt, region=(idx_0, idx_1))
            self.cxxclass.clear_buffer(self)

if pmi.isController:
    class DumpTopology(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.io.DumpTopologyLocal',
            pmicall=['dump', 'clear_buffer', 'observe_tuple', 'update'],
            pmiproperty=[],
            pmiinvoke=['get_data']
        )
