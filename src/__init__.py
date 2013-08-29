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
from espresso.FixedPairList import *
from espresso.FixedPairDistList import *
from espresso.FixedPairListAdress import *
from espresso.FixedTripleList import *
from espresso.FixedTripleAngleList import *
from espresso.FixedTripleListAdress import *
from espresso.FixedQuadrupleList import *
from espresso.FixedQuadrupleAngleList import *
from espresso.FixedTupleList import *
from espresso.FixedTupleListAdress import *
from espresso.Settle import *
from espresso.MultiSystem import *
from espresso.ParallelTempering import *
from espresso.Version import *
from espresso.CellList import *
from espresso.VirtualVerletList import *


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
