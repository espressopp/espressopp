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


# set up espresso basics
from espresso.main._setup import *

# load espresso into PMI
pmiimport('espresso')

import _espresso
from espresso.Exceptions import *
from espresso.Real3D import *
from espresso.RealND import *

from espresso.Tensor import *
from espresso.Int3D import *
from espresso.Particle import *
from espresso.ParticleGroup import *
from espresso.System import *
from espresso.VerletList import *
from espresso.VerletListTriple import *
from espresso.VerletListAdress import *
from espresso.FixedSingleList import *
from espresso.FixedPairList import *
from espresso.FixedPairDistList import *
from espresso.FixedPairListAdress import *
from espresso.FixedTripleList import *
from espresso.FixedTripleAngleList import *
from espresso.FixedTripleListAdress import *
from espresso.FixedQuadrupleList import *
from espresso.FixedQuadrupleListAdress import *
from espresso.FixedQuadrupleAngleList import *
from espresso.FixedTupleList import *
from espresso.FixedTupleListAdress import *
from espresso.MultiSystem import *
from espresso.ParallelTempering import *
from espresso.Version import *
from espresso.PLogger import *


infinity=float("inf")
nan=float("nan")
auto='auto'

# fetch the different subpackages
from espresso import esutil, bc, storage, integrator, interaction, analysis, tools, standard_system, external, check, io

if pmi.isController :
    # make sure that the workers exit when the script ends
    pmi.registerAtExit()
    # the script continues after this call
else :
    pmi.startWorkerLoop()
    # the script will usually not reach this point on the workers
