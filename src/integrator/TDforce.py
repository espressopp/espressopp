#  Copyright (C) 2017,2018
#      Jakub Krajniak (jkrajniak at gmail.com), Max Planck Institute for Polymer Research
#  Copyright (C) 2012,2013,2014,2015,2016
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
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
*****************************
espressopp.integrator.TDforce
*****************************

Thermodynamic force.

Example - how to turn on thermodynamic force (except for multiple moving spherical regions)

>>> fthd="tabletf.xvg"
>>> thdforce = espressopp.integrator.TDforce(system,verletlist) #info about centre and shape of adress region come from the verletlist, info about size of adress region not needed, tabulated file tabletf.xvg should be appropriate for the region size
>>> thdforce.addForce(itype=3,filename="tabletf.xvg",type=typeCG)
>>> integrator.addExtension(thdforce)

Example - how to turn on thermodynamic force for multiple moving spherical regions

>>> fthd="tabletf.xvg"
>>> thdforce = espressopp.integrator.TDforce(system, verletlist, startdist = 0.9, enddist = 2.1, edgeweightmultiplier = 20) #info about moving centres come from the verletlist. Info about size of adress region not needed, tabulated file tabletf.xvg should be appropriate for the region size (enddist - startdist)
>>> thdforce.addForce(itype=3,filename="tabletf.xvg",type=typeCG)
>>> integrator.addExtension(thdforce)

.. function:: espressopp.integrator.TDforce(system, verletlist, startdist, enddist, edgeweightmultiplier, slow)

        :param system: system object
        :param verletlist: verletlist object
        :param startdist: (default: 0.0) starting distance from center at which the TD force is actually applied. Needs to be altered when using several moving spherical regions (not used for static or single moving region)
        :param enddist: (default: 0.0) end distance from center up to which the TD force is actually applied. Needs to be altered when using several moving spherical regions (not used for static or single moving region)
        :param edgeweightmultiplier: (default: 20) interpolation parameter for multiple overlapping spherical regions (see Kreis et al., JCTC doi: 10.1021/acs.jctc.6b00440), the default should be fine for most applications (not used for static or single moving region)
        :param slow: (default: False) When used with RESPA Velocity Verlet, this flag decides whether the TD force is applied together with the slow, less frequently updated forces (slow=True) or with the fast, more frequently updated (slow=False) forces.
        :type system: shared_ptr<System>
        :type verletlist: shared_ptr<VerletListAdress>
        :type startdist: real
        :type enddist: real
        :type edgeweightmultiplier: int
        :type slow: bool

.. function:: espressopp.integrator.TDforce.addForce(itype, filename, type)

        Adds a thermodynamic force acting on particles of type "type".

        :param itype: interpolation type 1: linear, 2: Akima, 3: Cubic
        :param filename: filename for TD force file
        :param type: particle type on which the TD force needs to be applied
        :type itype: int
        :type filename: string
        :type type: int

.. function:: espressopp.integrator.TDforce.computeTDEnergy()

        Computes the energy corresponding to the thermodynamics force (summing over all different particle types and thermodynamic forces on them).

        :rtype: real
"""
from espressopp.esutil import cxxinit
from espressopp import pmi

from _espressopp import integrator_TDforce

class TDforceLocal(integrator_TDforce):

    def __init__(self, system, verletlist, startdist=0.0, enddist=0.0, edgeweightmultiplier=20, slow=False):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, integrator_TDforce, system, verletlist, startdist, enddist, edgeweightmultiplier, slow)

    def addForce(self, itype, filename, type):
            """
            Each processor takes the broadcasted interpolation type,
            filename and particle type
            """
            if pmi.workerIsActive():
                self.cxxclass.addForce(self, itype, filename, type)

    def computeTDEnergy(self):
            if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
              return self.cxxclass.computeTDEnergy(self)

if pmi.isController :
    class TDforce(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls =  'espressopp.integrator.TDforceLocal',
            pmiproperty = [ 'itype', 'filename'],
            pmicall = ['addForce', 'getForce', 'computeTDEnergy']
            )
