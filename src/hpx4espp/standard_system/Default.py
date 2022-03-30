#  Copyright (C) 2021-2022
#      Max Planck Institute for Polymer Research & JGU Mainz
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
**********************************
espressopp.hpx4espp.standard_system.Default
**********************************


.. py:method:: espressopp.hpx4espp.standard_system.Default(box, numSubs, rc = 1.12246, skin = 0.3, dt = 0.005, temperature = None)

                :param box:
                :param numSubs:
                :param real rc:
                :param real skin:
                :param real dt:
                :param temperature:
                :type box:
                :type temperature:

                Return default system and integrator, no interactions, no particles are set
                if tempearture is != None then Langevin thermostat is set to temperature (gamma is 1.0)
"""
import espressopp
import mpi4py.MPI as MPI

def Default(box, numSubs, numCommSubs=1, rc=1.12246, skin=0.3, dt=0.005, temperature=None, halfCellInt=1):

    assert(espressopp.hpx4espp.enabled())

    system      = espressopp.hpx4espp.System()
    system.rng  = espressopp.esutil.RNG()
    system.bc   = espressopp.bc.OrthorhombicBC(system.rng, box)
    system.skin = skin
    nodeGrid    = espressopp.tools.decomp.nodeGrid(MPI.COMM_WORLD.size,box,rc,skin)
    cellGrid    = espressopp.tools.decomp.cellGrid(box, nodeGrid, rc, skin, halfCellInt)

    # TODO: Set as default parameters in DomainDecomp
    commAsync = False
    excgAligned = False
    commUseChannels = False
    decompUseParFor = True

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--HPX4ESPP_OPTIMIZE_COMM", type=int, default=1)
    args, _ = parser.parse_known_args()
    HPX4ESPP_OPTIMIZE_COMM = args.HPX4ESPP_OPTIMIZE_COMM

    # setup DomainDecomp parameters
    subSize = [box[i]/nodeGrid[i] for i in range(3)]
    subNodeGrid = espressopp.hpx4espp.storage.nodeGridMultiple(numSubs, cellGrid)
    subCellGrid = espressopp.tools.decomp.cellGrid(subSize, subNodeGrid, rc, skin)

    expCellGrid = [subCellGrid[i]*subNodeGrid[i] for i in range(3)]
    for i in range(3):
        assert (expCellGrid[i] == cellGrid[i]), \
            "Mismatch in cellGrid. Expected {}, actual {}.\n" \
            "Make sure that subNodeGrid {} divides cellGrid {}. ".format(
                expCellGrid, cellGrid, subNodeGrid, cellGrid)

    system.storage = espressopp.hpx4espp.storage.DomainDecomposition(
        system, nodeGrid, cellGrid, halfCellInt, subCellGrid, numCommSubs, commAsync, excgAligned,
        commUseChannels, decompUseParFor, HPX4ESPP_OPTIMIZE_COMM)

    integrator     = espressopp.hpx4espp.integrator.VelocityVerlet(system)
    integrator.dt  = dt
    if (temperature != None):
        # determine expected number of threads from command line
        import argparse
        parser = argparse.ArgumentParser()
        parser.add_argument("--hpx:threads", type=int, dest="hpxThreads", required=True)
        args, _ = parser.parse_known_args()

        system.rngThread       = espressopp.hpx4espp.esutil.RNGThread(args.hpxThreads)
        thermostat             = espressopp.hpx4espp.integrator.LangevinThermostat(system)
        thermostat.gamma       = 1.0
        thermostat.temperature = temperature
        integrator.addExtension(thermostat)

    return system, integrator
