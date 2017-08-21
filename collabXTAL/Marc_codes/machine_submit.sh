#!/bin/bash

source /home/theorie/radumarc/Software/Espressopp/espressopp-1.9.2/ESPRC

ulimit -s unlimited

#mpirun -np 2 python LJ_PARRINELLO.py | tee logfile
#mpirun -np 2 python LJ_POLYMER.py | tee logfile
mpirun -np 2 python LJ_CMUMD.py | tee logfile
#mpirun -np 2 python LJ_CMUMD_POLYMER.py | tee logfile


