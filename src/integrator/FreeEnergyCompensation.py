#  Copyright (C) 2012-2018
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008-2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
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
********************************************
espressopp.integrator.FreeEnergyCompensation
********************************************

Free Energy Compensation used in Hamiltonian Adaptive Resolution Simulations (H-AdResS) or Path Integral Adaptive Resolution Simulations (PI-AdResS). This works for spherical or slab adaptive resolution geometries. However, it only works for fixed, non-moving atomistic region (otherwise, H-AdResS is not properly defined).

Example:

>>> fec = espressopp.integrator.FreeEnergyCompensation(system, center=[Lx/2, Ly/2, Lz/2])
>>> # set up the fec module with the center in the center of the box
>>> fec.addForce(itype=3,filename="tablefec.xvg",type=typeCG)
>>> # set up the actual force
>>> integrator.addExtension(fec)
>>> # add to previously defined integrator

.. function:: espressopp.integrator.FreeEnergyCompensation(system, center, sphereAdr, ntrotter, slow)

        :param system: system object
        :param center: (default: [], corresponds to (0.0, 0.0, 0.0) position) center of high resolution region
        :param sphereAdr: (default: False) Spherical AdResS region (True) vs. slab geometry with resolution change in x-direction (False)
        :param ntrotter: (default: 1) Trotter number when used in Path Integral AdResS. Default leads to normal non-PI-AdResS behaviour.
        :param slow: (default: False) When used with RESPA Velocity Verlet, this flag decides whether the Free Energy Compensation is applied together with the slow, less frequently updated forces (slow=True) or with the fast, more frequently updated (slow=False) forces.
        :type system: shared_ptr<System>
        :type center: list of reals
        :type sphereAdr: bool
        :type ntrotter: int
        :type slow: bool

.. function:: espressopp.integrator.FreeEnergyCompensation.addForce(itype, filename, type)

        :param itype: interpolation type 1: linear, 2: Akima, 3: Cubic
        :param filename: filename for TD force file
        :param type: particle type on which the TD force needs to be applied
        :type itype: int
        :type filename: string
        :type type: int

.. function:: espressopp.integrator.FreeEnergyCompensation.computeCompEnergy()

        :rtype: real
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_FreeEnergyCompensation

class FreeEnergyCompensationLocal(ExtensionLocal, integrator_FreeEnergyCompensation):

    def __init__(self, system, center=[], sphereAdr = False, ntrotter=1, slow=False):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_FreeEnergyCompensation, system, sphereAdr, ntrotter, slow)

            # set center of FreeEnergyCompensation force
            if (center != []):
                self.cxxclass.setCenter(self, center[0], center[1], center[2])

    def addForce(self, itype, filename, type):
            """
            Each processor takes the broadcasted interpolation type,
            filename and particle type
            """
            if pmi.workerIsActive():
                self.cxxclass.addForce(self, itype, filename, type)

    def computeCompEnergy(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.computeCompEnergy(self)

if pmi.isController :
    class FreeEnergyCompensation(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.FreeEnergyCompensationLocal',
            pmiproperty = [ 'itype', 'filename'],
            pmicall = ['addForce' , 'computeCompEnergy']
            )
