# Copyright (c) 2015-2018
#     Jakub Krajniak (jkrajniak at gmail.com)
#
# Copyright (c) 2015
#     Pierre de Buyl
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  ESPResSo++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


"""
**************************
**DumpH5MD** - IO object
**************************

This module provides a writer for H5MD_ file format.

.. function:: espressopp.io.DumpH5MD(system, integrator, filename, group_name, *args)

    :param system: The system object.
    :type system: espressopp.System
    :param integrator: System integrator.
    :type integrator: espressopp.integrator.MDIntegrator
    :param filename: The file name.
    :type filename: str
    :param group_name: The name of particle group.
    :type group_name: str
    :param is_adress: If positive then store position of AdResS particles
    :type is_adress: bool
    :param author: The name of author of this file
    :type author: str
    :param email: The e-mail address to the author.
    :type email: str
    :param chunk_size:
    :type chunk_size: int
    :param static_box: box size written as time-independent variable
    :type static_box: bool
    :param is_single_prec: Store float values with single precision (default: False)
    :type is_single_prec: bool
    :param store_position: Saves postions of particles
    :param store_species: Saves types of particles.
    :param store_state: Saves states of particles.
    :param store_velocity: Saves velocities of particles.
    :param store_force: Saves forces of particles
    :param store_charge: Saves charges of particles
    :param store_lambda: Saves lambdas (AdResS) of particles
    :param store_res_id: Saves residues id of particles.
    :param store_mass: Saves masses of particles

    :rtype: The DumpH5MD writer.

Example
+++++++

>>> traj_file = espressopp.io.DumpH5MD(
        system, integrator, output_file,
        group_name='atoms',
        static_box=False,
        author='xxx',
        email='xxx@xxx',
        store_species=True,
        store_velocity=True,
        store_state=True,
        store_lambda=True)

>>> for s in range(steps):
        integrator.run(int_steps)
        traj_file.dump(s*int_steps, s*int_steps*integrator.dt)

Important note. Within the current approach, this extension is not compatible with ExtAnalyze module.
Therefore, this code does not work:

>>> ext_analyze = espressopp.integrator.ExtAnalyze(traj_file, 10)
>>> integrator.addExtension(ext_analyze)
>>> integrator.run(2000)

.. _sorting-file-label:

Sorting file
++++++++++++++

The content of the `/particles/{}/` is not sorted with respect to the
particle id. This is because of the way how the data are stored
by multiple cores simultaneously.

If the flag `do_sort` is True then during the close method, the data
will be sorted.


.. _H5MD: http://nongnu.org/h5md/

"""

import espressopp
from espressopp.esutil import cxxinit
from espressopp import pmi
from _espressopp import io_DumpH5MD
from mpi4py import MPI
import h5py
import numpy as np
import pyh5md
import os

import time as py_time


class DumpH5MDLocal(io_DumpH5MD):
    def __init__(self, system, integrator, filename, group_name='atoms',
                 store_position=True,
                 store_species=True,
                 store_state=False,
                 store_velocity=False,
                 store_force=False,
                 store_charge=False,
                 store_lambda=False,
                 store_res_id=False,
                 store_mass=True,
                 static_box=True,
                 is_adress=False,
                 author='xxx',
                 email='xxx',
                 chunk_size=256,
                 is_single_prec=False,
                 do_sort=True):
        """
        Args:
            system: The system object.
            integrator: The integrator.
            filename: The name of hdf file name.
            group_name: The name of atom groups. (default: 'atoms').
            store_position: If set to True then position will be stored. (default: True)
            store_species: If set to True then species will be stored. (default: True)
            store_state: If set to True then state will be stored. (default: False)
            store_velocity: If set to True then velocity will be stored. (default: False)
            store_force: If set to True then force will be stored. (default: False)
            store_charge: If set to True then charge will be stored. (default: False)
            store_lambda: If set to True then lambda (AdResS) will be stored. (default: False)
            store_res_id: If set to True then store res_id. (default: False)
            store_mass: If set to True then store mass. (default: True)
            static_box: If set to True then box is static (like in NVT ensemble) (default: True)
            is_adress: If set to True then AdResS particles will be save instead of
                coarse-grained.
            author: The name of author of the file. (default: xxx)
            email: The e-mail to author of that file. (default: xxx)
            chunk_size: The size of data chunk. (default: 128)
            is_single_prec: Use single precision instead of double.
            do_sort: If set to True then HDF5 will be sorted on close.
        """
        if not pmi.workerIsActive():
            return
        cxxinit(self, io_DumpH5MD, system, is_adress)

        self.filename = filename
        self.group_name = group_name
        self.store_position = store_position
        self.store_species = store_species
        self.store_state = store_state
        self.store_velocity = store_velocity
        self.store_force = store_force
        self.store_charge = store_charge
        self.store_lambda = store_lambda
        self.store_res_id = store_res_id
        self.store_mass = store_mass
        self.static_box = static_box
        self.chunk_size = chunk_size
        self.do_sort = do_sort
        self.single_prec = is_single_prec

        self.system = system
        self.integrator = integrator

        if pmi.isController:
            if os.path.exists(filename):
                basename = os.path.basename(filename)
                dirname = os.path.dirname(filename)
                new_filename = '{}/{}_{}'.format(dirname, int(py_time.time()), os.path.basename(filename))
                os.rename(filename, new_filename)
                print('File {} exists, moved to {}'.format(filename, new_filename))
        try: 
            self.file = pyh5md.File(
                filename, 'w',
                creator='espressopp',
                creator_version=espressopp.VersionLocal().info(),
                author=author,
                author_email=email,
                driver='mpio',
                comm=MPI.COMM_WORLD)
        except NameError:
            # h5py in serial mode.
            self.file = pyh5md.File(
                filename, 'w',
                creator='espressopp',
                creator_version=espressopp.VersionLocal().info(),
                author=author,
                author_email=email)

        self._system_data()

        self.float_type = np.float32 if is_single_prec else np.float64
        self.int_type = np.int32 if is_single_prec else np.int

        part = self.file.particles_group(self.group_name)
        self.particle_group = part

        if self.static_box:
            part.create_box(
                dimension=3,
                boundary=['periodic', 'periodic', 'periodic'],
                store='fixed',
                data=np.array([ed_i for ed_i in self.system.bc.boxL], dtype=self.float_type))
        else:
            part.create_box(
                dimension=3,
                boundary=['periodic', 'periodic', 'periodic'],
                store='time', time=True,
                data=np.zeros(3, dtype=self.float_type))

        self.id_e = pyh5md.element(part, 'id', store='time', time=True, shape=(self.chunk_size,), maxshape=(None, ),
                                   dtype=self.int_type, fillvalue=-1)
        if store_mass:
            self.mass = pyh5md.element(part, 'mass', store='time', time=True, maxshape=(None, ),
                                       shape=(self.chunk_size,), dtype=self.float_type, fillvalue=-1)
        if self.store_position:
            self.position = pyh5md.element(part, 'position', store='time', time=True, maxshape=(None, 3),
                                           shape=(self.chunk_size, 3), dtype=self.float_type)
            self.image = pyh5md.element(part, 'image', store='time', time=True, maxshape=(None, 3),
                                        shape=(self.chunk_size, 3), dtype=self.float_type)
        if self.store_species:
            self.species = pyh5md.element(part, 'species', store='time', time=True, maxshape=(None, ),
                                          shape=(self.chunk_size,), dtype=self.int_type, fillvalue=-1)
        if self.store_state:
            self.state = pyh5md.element(part, 'state', store='time', time=True, maxshape=(None, ),
                                        shape=(self.chunk_size,), dtype=self.int_type,  fillvalue=-1)
        if self.store_velocity:
            self.velocity = pyh5md.element(part, 'velocity', store='time', time=True, maxshape=(None, 3),
                                           shape=(self.chunk_size, 3), dtype=self.float_type)
        if self.store_force:
            self.force = pyh5md.element(part, 'force', store='time', time=True, maxshape=(None, 3),
                                        shape=(self.chunk_size, 3), dtype=self.float_type)
        if self.store_charge:
            self.charge = pyh5md.element(part, 'charge', store='time', time=True, maxshape=(None, ),
                                         shape=(self.chunk_size,), dtype=self.float_type, fillvalue=-1)
        if self.store_lambda:
            self.lambda_adr = pyh5md.element(part, 'lambda_adr', store='time', time=True, maxshape=(None, ),
                                             shape=(self.chunk_size,), dtype=self.float_type,  fillvalue=-1)
        if self.store_res_id:
            self.res_id = pyh5md.element(part, 'res_id', store='time', time=True, maxshape=(None, ),
                                         shape=(self.chunk_size, ), dtype=self.int_type,  fillvalue=-1)

        self.commTimer = 0.0
        self.updateTimer = 0.0
        self.writeTimer = 0.0
        self.flushTimer = 0.0
        self.closeTimer = 0.0
        self.resizeCounter = 0

    def _system_data(self):
        """Stores specific information about simulation."""
        # Creates /system group
        parameters = {
            'software-id': 'espressopp',
            'rng-seed': self.system.rng.get_seed(),
            'skin': self.system.skin,
        }
        parameters['dt'] = self.integrator.dt
        self.set_parameters(parameters)

    def getTimers(self):
        if pmi.workerIsActive():
            return {'commTimer': self.commTimer,
                    'updateTimer': self.updateTimer,
                    'writeTimer': self.writeTimer,
                    'flushTimer': self.flushTimer,
                    'closeTimer': self.closeTimer,
                    'resizeCounter': self.resizeCounter
                   }

    def set_parameters(self, paramters):
        if pmi.workerIsActive():
            if 'parameters' not in self.file:
                self.file.create_group('parameters')
            g_params = self.file['parameters']
            for k, v in paramters.iteritems():
                g_params.attrs[k] = v

    def get_file(self):
        if pmi.workerIsActive():
            return self.file

    def update(self):
        if pmi.workerIsActive():
            self.cxxclass.update(self)

    def clear_buffers(self):
        if pmi.workerIsActive():
            self.cxxclass.clear_buffers(self)

    def getPosition(self):
        if pmi.workerIsActive():
            return self.cxxclass.getPosition(self)

    def getImage(self):
        if pmi.workerIsActive():
            return self.cxxclass.getImage(self)

    def getVelocity(self):
        if pmi.workerIsActive():
            return self.cxxclass.getVelocity(self)

    def getForce(self):
        if pmi.workerIsActive():
            return self.cxxclass.getForce(self)

    def getId(self):
        if pmi.workerIsActive():
            return self.cxxclass.getId(self)

    def getSpecies(self):
        if pmi.workerIsActive():
            return self.cxxclass.getSpecies(self)

    def getState(self):
        if pmi.workerIsActive():
            return self.cxxclass.getState(self)

    def getCharge(self):
        if pmi.workerIsActive():
            return self.cxxclass.getCharge(self)

    def getMass(self):
        if pmi.workerIsActive():
            return self.cxxclass.getMass(self)

    def getLambdaAdr(self):
        if pmi.workerIsActive():
            return self.cxxclass.getLambda(self)

    def getResId(self):
        if pmi.workerIsActive():
            return self.cxxclass.getResId(self)

    def dump(self, step=0, time=0):
        if not pmi.workerIsActive():
            return

        # Take it directly from integrator;
        step = self.integrator.step
        time = step*self.integrator.dt

        time0 = py_time.time()
        self.update()
        self.updateTimer += (py_time.time() - time0)

        time0 = py_time.time()
        NLocal = np.array(self.NLocal, 'i')
        NMaxLocal = np.array(0, 'i')
        MPI.COMM_WORLD.Allreduce(NLocal, NMaxLocal, op=MPI.MAX)
        cpu_size = ((NMaxLocal//self.chunk_size)+1)*self.chunk_size
        total_size = MPI.COMM_WORLD.size*cpu_size
        idx_0 = MPI.COMM_WORLD.rank*cpu_size
        idx_1 = idx_0+NLocal
        self.commTimer += (py_time.time() - time0)

        collective_write = False
        isResized = False

        time0 = py_time.time()
        # Store ids. Always!
        id_ar = np.asarray(self.getId())
        if total_size > self.id_e.value.shape[1]:
            isResized = True
            self.id_e.value.resize(total_size, axis=1)
        self.id_e.append(id_ar, step, time, region=(idx_0, idx_1), collective=collective_write)

        # Store box values at every time step
        if not self.static_box:
            self.particle_group.box.edges.append(
                np.array([edge_i for edge_i in self.system.bc.boxL], dtype=self.float_type),
                step,
                time, collective=collective_write)

        if self.store_position:
            pos_data = self.getPosition()
            pos = np.asarray(pos_data, dtype=self.float_type)
            if total_size > self.position.value.shape[1]:
                isResized = True
                self.position.value.resize(total_size, axis=1)
            self.position.append(pos, step, time, region=(idx_0, idx_1), collective=collective_write)
            # Store image.
            data_image = self.getImage()
            if data_image is None:
                image = np.asarray(None, dtype=self.float_type)
            else:
                image = np.asarray(data_image, dtype=self.int_type)
            if total_size > self.image.value.shape[1]:
                isResized = True
                self.image.value.resize(total_size, axis=1)
            self.image.append(image, step, time, region=(idx_0, idx_1), collective=collective_write)

        # Store velocity.
        if self.store_velocity:
            vel = np.asarray(self.getVelocity(), dtype=self.float_type)
            if total_size > self.velocity.value.shape[1]:
                isResized = True
                self.velocity.value.resize(total_size, axis=1)
            self.velocity.append(vel, step, time, region=(idx_0, idx_1), collective=collective_write)

        if self.store_force:
            force = np.asarray(self.getForce(), dtype=self.float_type)
            if total_size > self.force.value.shape[1]:
                isResized = True
                self.force.value.resize(total_size, axis=1)
            self.force.append(force, step, time, region=(idx_0, idx_1), collective=collective_write)

        if self.store_charge:
            charge = np.asarray(self.getCharge(), dtype=self.float_type)
            if total_size > self.charge.value.shape[1]:
                isResized = True
                self.charge.value.resize(total_size, axis=1)
            self.charge.append(charge, step, time, region=(idx_0, idx_1), collective=collective_write)

        # Store mass.
        if self.store_mass:
            mass = np.asarray(self.getMass(), dtype=self.float_type)
            if total_size > self.mass.value.shape[1]:
                isResized = True
                self.mass.value.resize(total_size, axis=1)
            self.mass.append(mass, step, time, region=(idx_0, idx_1), collective=collective_write)

        # Store species.
        if self.store_species:
            species = np.asarray(self.getSpecies())
            if total_size > self.species.value.shape[1]:
                isResized = True
                self.species.value.resize(total_size, axis=1)
            self.species.append(species, step, time, region=(idx_0, idx_1), collective=collective_write)

        # Store state.
        if self.store_state:
            state = np.asarray(self.getState())
            if total_size > self.state.value.shape[1]:
                isResized = True
                self.state.value.resize(total_size, axis=1)
            self.state.append(state, step, time, region=(idx_0, idx_1), collective=collective_write)

        # Store lambda_adr
        if self.store_lambda:
            lambda_adr = np.asarray(self.getLambda())
            if total_size > self.lambda_adr.value.shape[1]:
                isResized = True
                self.lambda_adr.value.resize(total_size, axis=1)
            self.lambda_adr.append(lambda_adr, step, time, region=(idx_0, idx_1), collective=collective_write)

        # Store res_id
        if self.store_res_id:
            res_id = np.asarray(self.getResId())
            if total_size > self.res_id.value.shape[1]:
                isResized = True
                self.res_id.value.resize(total_size, axis=1)
            self.res_id.append(res_id, step, time, region=(idx_0, idx_1), collective=collective_write)
        self.writeTimer += (py_time.time() - time0)

        if isResized:
            self.resizeCounter += 1

    def close(self):
        if pmi.workerIsActive():
            time0 = py_time.time()
            self.file.close()
            self.closeTimer += (py_time.time() - time0)

    def flush(self):
        if pmi.workerIsActive():
            time0 = py_time.time()
            self.file.flush()
            self.flushTimer += (py_time.time() - time0)


if pmi.isController:
    def sort_file(h5):
        """Sort data file."""
        atom_groups = [ag for ag in h5['/particles'] if 'id' in h5['/particles/{}/'.format(ag)]]
        T = len(h5['/particles/{}/id/value'.format(atom_groups[0])])
        # Iterate over time frames.
        for t in xrange(T):
            for ag in atom_groups:
                ids = h5['/particles/{}/id/value'.format(ag)]
                idd = [
                    x[1] for x in sorted(
                        [(p_id, col_id) for col_id, p_id in enumerate(ids[t])],
                        key=lambda y: (True, y[0]) if y[0] == -1 else (False, y[0]))
                    ]
                for k in h5['/particles/{}/'.format(ag)].keys():
                    if 'value' in h5['/particles/{}/{}'.format(ag, k)].keys():
                        path = '/particles/{}/{}/value'.format(ag, k)
                        h5[path][t] = h5[path][t][idd]

    class DumpH5MD(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls='espressopp.io.DumpH5MDLocal',
            pmicall=['update', 'getPosition', 'getId', 'getSpecies', 'getState', 'getImage',
                     'getVelocity', 'getMass', 'getCharge', 'getResId',
                     'dump', 'clear_buffers', 'flush', 'get_file', 'set_parameters'],
            pmiinvoke=['getTimers'],
            pmiproperty=['store_position', 'store_species', 'store_state', 'store_velocity',
                         'store_charge', 'store_res_id', 'store_lambda', 'store_mass'])

        def close(self):
            pmi.call(self.pmiobject, "close")
            # Sort file if flag is set to true.
            if self.pmiobject.do_sort:
                h5 = h5py.File(self.pmiobject.filename, 'r+')
                print('Sorting file, please wait...')
                sort_file(h5)
                print('File sorted')
                h5.close()
