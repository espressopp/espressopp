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


r"""
**********************
espressopp.esutil.Grid
**********************


"""
from espressopp import pmi
from _espressopp import esutil_Grid

class GridLocal(esutil_Grid):
  pass

if pmi.isController:
    class Grid(object, metaclass=pmi.Proxy):
        'Grid class'
        pmiproxydefs = dict(
            cls = 'espressopp.esutil.GridLocal',
            localcall = [ 'mapIndexToPosition' ]
            #localcall = [ '__call__', 'normal', 'gamma', 'uniformOnSphere' ],
            #pmicall = [ 'seed' ]
        )
    
