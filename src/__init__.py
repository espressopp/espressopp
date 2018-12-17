#  Copyright (C) 2012,2013
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


# set up espressopp basics
from espressopp.main._setup import *

# load espressopp into PMI
pmiimport('espressopp')

import _espressopp
from espressopp.Exceptions import *
from espressopp.Real3D import *
from espressopp.Quaternion import *
from espressopp.RealND import *

from espressopp.Tensor import *
from espressopp.Int3D import *
from espressopp.Particle import *
from espressopp.ParticleGroup import *
from espressopp.System import *
from espressopp.VerletList import *
from espressopp.VerletListTriple import *
from espressopp.VerletListAdress import *
from espressopp.FixedSingleList import *
from espressopp.FixedPairList import *
from espressopp.FixedPairDistList import *
from espressopp.FixedPairListAdress import *
from espressopp.FixedTripleList import *
from espressopp.FixedTripleAngleList import *
from espressopp.FixedTripleListAdress import *
from espressopp.FixedQuadrupleList import *
from espressopp.FixedQuadrupleListAdress import *
from espressopp.FixedQuadrupleAngleList import *
from espressopp.FixedTupleList import *
from espressopp.FixedTupleListAdress import *
from espressopp.FixedLocalTupleList import *
from espressopp.MultiSystem import *
from espressopp.ParallelTempering import *
from espressopp.Version import *
from espressopp.PLogger import *


infinity=float("inf")
nan=float("nan")
auto='auto'

# fetch the different subpackages
from espressopp import esutil, bc, storage, integrator, interaction, analysis, tools, standard_system, external, check, io

if pmi.isController :
    # make sure that the workers exit when the script ends
    pmi.registerAtExit()
    # the script continues after this call
else :
    pmi.startWorkerLoop()
    # the script will usually not reach this point on the workers
