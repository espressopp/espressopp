# set up espresso paths and logging
from espresso.main import __setup

from espresso.esutil import pmiimport
pmiimport('espresso')

from espresso.Real3D import *
from espresso.Int3D import *
from espresso import pmi

infinity=float("inf")
nan=float("nan")

if pmi.isController :
    # make sure that the workers exit when the script ends
    pmi.registerAtExit()
    # the script continues after this call
else :
    pmi.startWorkerLoop()
    # the script will usually not reach this point on the workers


    
