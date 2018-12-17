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
****************************
espressopp.integrator.Adress
****************************

The AdResS object is an extension to the integrator. It makes sure that the
integrator also processes the atomistic particles and not only the CG particles.
Hence, this object is of course only used when performing AdResS or H-AdResS
simulations.

In detail the AdResS extension makes sure:

* that also the forces on the atomistic particles are initialized and set to
  by Adress::initForces
* that also the atomistic particles are integrated and propagated by
  Adress::integrate1 and Adress::integrate2

Example - how to turn on the AdResS integrator extension:

>>> adress = espressopp.integrator.Adress(system, verletlist, fixedtuplelist)
>>> integrator.addExtension(adress)

If KTI is set to True, then the resolution parameters are not updated. This can be used for example for Kirkwood thermodynamic integration, during which one manually sets the whole system on different resolution parameters. KTI = True then prevents overwriting these manually set values.
Furthermore, when having moving AdResS regions based on particles, regionupdates specifies the update frequency of the AdResS region in number of steps (or, to be more precise, calls of communicateAdrPositions()). Note that there is a tradeoff: The more frequently the AdResS region is updated, the more gradually and accurately the AdResS region changes and adapts it shape. This could allow for a smaller overall AdResS region and possibly a smoother simulation. However, when having many AdResS region defining particles, these frequent updates can become computationally significant and cost additional simulation time. The optimum is highly system and application dependent.

Finally, when making use of a RESPA Velocity Verlet integrator, then the multistep parameter defines after how many steps of the inner integration loop the slow forces are updated. It should be set consistently with the same parameter in VelocityVerletRESPA.

.. py:class:: espressopp.integrator.Adress(_system, _verletlist, _fixedtuplelist, KTI, regionupdates, multistep)

                :param _system: system object
                :param _verletlist: verletlist object
                :param _fixedtuplelist: fixedtuplelist object
                :param KTI: (default: False) update resolution parameter? (Yes: set False, No: set True)
                :param regionupdates: (default: 1) after how many steps does the AdResS region needs to be updated?
                :param multistep: (default: 1) when used with VelocityVerletRESPA (otherwise, ignored), after how many steps of the inner integration loop do we update the slow forces? This parameter should be set consistently with multistep in VelocityVerletRESPA.
                :type _system: shared_ptr<System>
                :type _verletlist: shared_ptr<VerletListAdress>
                :type _fixedtuplelist: shared_ptr<FixedTupleListAdress>
                :type KTI: bool
                :type regionupdates: int
                :type multistep: int
"""

from espressopp.esutil import cxxinit
from espressopp import pmi

from espressopp.integrator.Extension import *
from _espressopp import integrator_Adress

class AdressLocal(ExtensionLocal, integrator_Adress):


    def __init__(self, _system, _verletlist, _fixedtuplelist, KTI = False, regionupdates = 1, multistep = 1):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            if pmi.workerIsActive():
                cxxinit(self, integrator_Adress, _system, _verletlist, _fixedtuplelist, KTI, regionupdates, multistep)

if pmi.isController:
    class Adress(Extension):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espressopp.integrator.AdressLocal' #,
            #pmiproperty = [ 'builds' ],
            #pmicall = [ 'totalSize', 'exclude', 'addAdrParticles', 'rebuild' ]
            )
