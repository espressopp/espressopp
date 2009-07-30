# set up espresso paths and logging
from espresso import __setup

from espresso.esutil import pmiimport
pmiimport('espresso')

from espresso.Property import *
from espresso.Real3D import *
from espresso import pmi

if pmi.IS_CONTROLLER :
    # make sure that the workers exit when the script ends
    pmi.registerAtExit()
    # the script continues after this call
else :
    pmi.startWorkerLoop()
    # the script will usually not reach this point on the workers


    
