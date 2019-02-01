#  Copyright (c) 2015-2018
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
**DumpTopology** - IO object
*********************************************

.. function:: espressopp.io.DumpTopology(system, integrator, h5md_file)

   :param system: The ESPP System object.
   :type system: espressopp.System
   :param integrator: The integrator object.
   :type integrator: espressopp.integrator.MDIntegrator
   :param h5md_file: The H5MD file object.
   :type h5md_file: espressopp.io.DumpH5MD

.. function:: espressopp.io.DumpTopology.dump()

   Store data from tuple into memory.

.. function:: espressopp.io.DumpTopology.observe_tuple(fpl, name, particle_group)

   :param fpl: The FixedPairList object.
   :type fpl: espressopp.FixedPairList
   :param name: The name of the tuple to store in H5MD file.
   :type name: str
   :param particle_group: The particle group to referee to.
   :type particle_group: str

.. function:: espressopp.io.DumpTopology.update()

   Update the H5MD file.

.. function:: espressopp.io.DumpTopology.add_static_tuple(fpl, name, particle_group)

   Write data from fixed pair list once, not updated whenever the tuple is changed.

   :param fpl: The FixedPairList object.
   :type fpl: espressopp.FixedPairList
   :param name: The name of the tuple to store in H5MD file.
   :type name: str
   :param particle_group: The particle group to referee to.
   :type particle_group: str

Example
+++++++

The code belows dump topology every 10 time steps and stores pairs from
FixedPairList `fpl`.

>>> traj_file = espressopp.io.DumpH5MD(
        system, output_file,
        group_name='atoms',
        static_box=False,
        author='xxx',
        email='xxx@xxx',
        store_species=True,
        store_velocity=True,
        store_state=True,
        store_lambda=True)

>>> dump_topol = espressopp.io.DumpTopology(system, integrator, traj_file)
>>> dump_topol.observe_tuple(fpl_a_a, 'fpl')
>>> dump_topol.dump()
>>> dump_topol.update()
>>> ext_dump = espressopp.integrator.ExtAnalyze(dump_topol, 10)
>>> integrator.addExtension(ext_dump)

Stores static data from FixedPairList `fpl_0`, those are time-independent.

>>> dump_topol.add_static_tuple(fpl_0, 'fpl_0', 'atoms')

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
            self.triple_index = 0
            self.quadruple_index = 0
            self.tuple_data = {}
            self.triple_data = {}
            self.quadruple_data = {}
            if 'connectivity' not in self.h5md_file.file:
                self.h5md_file.file.create_group('connectivity')
            self.connectivity = self.h5md_file.file['connectivity']
            self.chunk_size = h5md_file.chunk_size
            self.dt = integrator.dt

    def dump(self):
        """Dump data to the internal buffer."""
        if pmi.workerIsActive():
            self.cxxclass.dump(self)

    def observe_tuple(self, fpl, name, particle_group='atoms'):
        if pmi.workerIsActive():
            self.cxxclass.observe_tuple(self, fpl)
            g = pyh5md.element(
                self.connectivity,
                name,
                store='time',
                maxshape=(None, 2),
                shape=(self.chunk_size, 2),
                dtype=self.h5md_file.int_type,
                fillvalue=-1)
            g.attrs['particle_group'] = particle_group
            self.tuple_data[self.tuple_index] = g
            self.tuple_index += 1

    def observe_triple(self, ftl, name, particle_group='atoms'):
        if pmi.workerIsActive():
            self.cxxclass.observe_triple(self, ftl)
            g = pyh5md.element(self.connectivity, name, store='time', maxshape=(None, 3), shape=(self.chunk_size, 3),
                               dtype=self.h5md_file.int_type, fillvalue=-1)
            g.attrs['particle_group'] = particle_group
            self.triple_data[self.triple_index] = g
            self.triple_index += 1

    def observe_quadruple(self, fql, name, particle_group='atoms'):
        if pmi.workerIsActive():
            self.cxxclass.observe_quadruple(self, fql)
            g = pyh5md.element(
                self.connectivity,
                name,
                store='time',
                maxshape=(None, 4),
                shape=(self.chunk_size, 4),
                dtype=self.h5md_file.int_type,
                fillvalue=-1)
            g.attrs['particle_group'] = particle_group
            self.quadruple_data[self.quadruple_index] = g
            self.quadruple_index += 1

    def add_static_tuple(self, fpl, name, particle_group='atoms'):
        if pmi.workerIsActive():
            bonds = fpl.getBonds()
            NMaxLocal = np.array(len(bonds), 'i')
            NMaxGlobal = np.array(0, 'i')
            MPI.COMM_WORLD.Allreduce(NMaxLocal, NMaxGlobal, op=MPI.MAX)
            size_per_cpu = ((NMaxGlobal // self.chunk_size)+1)*self.chunk_size
            total_size = MPI.COMM_WORLD.size*size_per_cpu
            # Prepares Dataset.
            g = pyh5md.element(
                self.connectivity,
                name,
                store='fixed',
                dtype=self.h5md_file.int_type,
                fillvalue=-1,
                shape=(total_size, 2))
            g.attrs['particle_group'] = particle_group
            # Calculates index per cpu.
            idx_0 = MPI.COMM_WORLD.rank*size_per_cpu
            idx_1 = idx_0 + NMaxLocal
            # Writes data.
            g[idx_0:idx_1] = bonds

    def add_static_triple(self, ftl, name, particle_group='atoms'):
        if pmi.workerIsActive():
            triplets = ftl.getTriples()
            NMaxLocal = np.array(len(triplets), 'i')
            NMaxGlobal = np.array(0, 'i')
            MPI.COMM_WORLD.Allreduce(NMaxLocal, NMaxGlobal, op=MPI.MAX)
            size_per_cpu = ((NMaxGlobal // self.chunk_size)+1)*self.chunk_size
            total_size = MPI.COMM_WORLD.size*size_per_cpu
            g = pyh5md.element(
                self.connectivity,
                name,
                store='fixed',
                dtype=self.h5md_file.int_type,
                fillvalue=-1,
                shape=(total_size, 3))
            g.attrs['particle_group'] = particle_group
            idx_0 = MPI.COMM_WORLD.rank*size_per_cpu
            idx_1 = idx_0 + NMaxLocal
            g[idx_0:idx_1] = triplets

    def add_static_quadruple(self, fql, name, particle_group='atoms'):
        if pmi.workerIsActive():
            quadruplets = fql.getQuadruples()
            NMaxLocal = np.array(len(quadruplets), 'i')
            NMaxGlobal = np.array(0, 'i')
            MPI.COMM_WORLD.Allreduce(NMaxLocal, NMaxGlobal, op=MPI.MAX)
            size_per_cpu = ((NMaxGlobal // self.chunk_size)+1)*self.chunk_size
            total_size = MPI.COMM_WORLD.size*size_per_cpu
            g = pyh5md.element(
                self.connectivity,
                name,
                store='fixed',
                dtype=self.h5md_file.int_type,
                fillvalue=-1,
                shape=(total_size, 4))
            g.attrs['particle_group'] = particle_group
            idx_0 = MPI.COMM_WORLD.rank*size_per_cpu
            idx_1 = idx_0 + NMaxLocal
            g[idx_0:idx_1] = quadruplets

    def update(self):
        """Load data from the buffer and store in the HDF5 file."""
        if pmi.workerIsActive():
            raw_data = self.cxxclass.get_data(self)
            step_data = {}
            max_sizes = [0, 0, 0] # {2: 0, 3: 0, 4: 0}
            while raw_data != []:
                step = raw_data.pop()
                fpl_idx = raw_data.pop()
                fpl_type = raw_data.pop()
                fpl_size = raw_data.pop()
                data = []
                for _ in range(fpl_size):
                    if fpl_type == 2:
                        b1 = raw_data.pop()
                        b2 = raw_data.pop()
                        data.append((b1, b2))
                    elif fpl_type == 3:
                        b1 = raw_data.pop()
                        b2 = raw_data.pop()
                        b3 = raw_data.pop()
                        data.append((b1, b2, b3))
                    elif fpl_type == 4:
                        b1 = raw_data.pop()
                        b2 = raw_data.pop()
                        b3 = raw_data.pop()
                        b4 = raw_data.pop()
                        data.append((b1, b2, b3, b4))
                    else:
                        raise RuntimeError('Wrong fpl_type: {}'.format(fpl_type))
                max_sizes[fpl_type-2] = max(len(data), max_sizes[fpl_type-2])
                if step not in step_data:
                    step_data[step] = {}
                if fpl_type not in step_data[step]:
                    step_data[step][fpl_type] = {}
                step_data[step][fpl_type][fpl_idx] = np.array(data, dtype=self.h5md_file.int_type)

            MPI.COMM_WORLD.Barrier()

            NMaxLocal = np.array(max_sizes, 'i')
            NMaxGlobal = np.zeros(3, 'i')
            MPI.COMM_WORLD.Allreduce(NMaxLocal, NMaxGlobal, op=MPI.MAX)
            cpu_size = ((NMaxGlobal // self.chunk_size)+1)*self.chunk_size
            total_size = MPI.COMM_WORLD.size*cpu_size
            idx_0 = MPI.COMM_WORLD.rank*cpu_size
            idx_1 = idx_0 + NMaxLocal
            for step in sorted(step_data):
                for fpl_type in step_data[step]:
                    for fpl_idx in step_data[step][fpl_type]:
                        data = step_data[step][fpl_type][fpl_idx]
                        if data.shape[0] < NMaxLocal[fpl_type-2]:
                            old_shape = data.shape
                            if len(data):
                                linear_data = data.reshape((old_shape[0]*old_shape[1],))
                            else:
                                linear_data = data
                            data = np.pad(linear_data, (0, fpl_type*NMaxLocal[fpl_type-2]-linear_data.shape[0]),
                                          'constant', constant_values=-1)
                            data = data.reshape((NMaxLocal[fpl_type-2], fpl_type))
                        if fpl_type == 2:
                            ds = self.tuple_data[fpl_idx]
                        elif fpl_type == 3:
                            ds = self.triple_data[fpl_idx]
                        elif fpl_type == 4:
                            ds = self.quadruple_data[fpl_idx]
                        else:
                            raise RuntimeError('Wrong fpl_type: {}'.format(fpl_type))
                        if total_size[fpl_type-2] > ds.value.shape[1]:
                            ds.value.resize(total_size[fpl_type-2], axis=1)
                        ds.append(data, step, step*self.dt, region=(idx_0[fpl_type-2], idx_1[fpl_type-2]))
            self.cxxclass.clear_buffer(self)

if pmi.isController:
    class DumpTopology(ParticleAccess):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.io.DumpTopologyLocal',
            pmicall=['dump', 'clear_buffer', 'observe_tuple', 'observe_triple', 'observe_quadruple', 'update',
                     'add_static_tuple', 'add_static_triple', 'add_static_quadruple'],
            pmiproperty=[],
            pmiinvoke=['get_data']
        )
