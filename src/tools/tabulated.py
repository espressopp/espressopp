#  Copyright (C) 2012,2013,2016
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

"""
********************************
tabulated - write tabulated file
********************************
"""

import espressopp
from espressopp import Real3D

def writeTabFile(pot, name, N, low=0.0, high=2.5, body=2):
    """
    writeTabFile can be used to create a table for any potential
    Parameters are:
    * pot     : this is any espressopp.interaction potential
    * name    : filename
    * N       : number of line to write
    * low     : lowest r (default is 0.0)
    * high    : highest r (default is 2.5)
    
    This function has not been tested for 3 and 4 body interactions
    """
    outfile = open(name, "w")
    delta = (high - low) / (N - 1)

    for i in xrange(N):
        r = low + i * delta
        energy = pot.computeEnergy(r)
        if body == 2:# this is for 2-body potentials
            force = pot.computeForce(Real3D(r, 0.0, 0.0))[0]
            #force /= r
        else: # this is for 3- and 4-body potentials
            force = pot.computeForce(r)
        outfile.write("%15.8g %15.8g %15.8g\n"%(r, energy, force))

    outfile.close()
