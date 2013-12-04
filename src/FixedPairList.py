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


"""
**************************
**espresso.FixedPairList**
**************************

"""
from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit
from math import sqrt

class FixedPairListLocal(_espresso.FixedPairList):
    'The (local) fixed pair list.'

    def __init__(self, storage):
        'Local construction of a fixed pair list'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.FixedPairList, storage)

    def add(self, pid1, pid2):
        'add pair to fixed pair list'
        if pmi.workerIsActive():
            return self.cxxclass.add(self, pid1, pid2)

    def size(self):
        'count number of bonds in GlobalPairList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.size(self)

    def addBonds(self, bondlist):
        """
        Each processor takes the broadcasted bondlist and
        adds those pairs whose first particle is owned by
        this processor.
        """
        
        if pmi.workerIsActive():
            for bond in bondlist:
                pid1, pid2 = bond
                self.cxxclass.add(self, pid1, pid2)

    def getBonds(self):
        'return the bonds of the GlobalPairList'
        if pmi.workerIsActive():
          bonds=self.cxxclass.getBonds(self)
          return bonds
      
    def resetLongtimeMaxBond(self):
        'reset long time maximum bond to 0.0' 
        if pmi.workerIsActive():
          self.cxxclass.resetLongtimeMaxBondSqr(self)
          
    def getLongtimeMaxBondLocal(self):
        'return the maximum bond length this pairlist ever had (since reset or construction)'
        if pmi.workerIsActive(): 
            mxsqr = self.cxxclass.getLongtimeMaxBondSqr(self)
            return sqrt(mxsqr)
            
if pmi.isController:
    class FixedPairList(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.FixedPairListLocal',
            #localcall = [ 'add' ],
            pmicall = [ 'add', 'addBonds', 'resetLongtimeMaxBond' ],
            pmiinvoke = ['getBonds', 'size', 'getLongtimeMaxBondLocal']
        )
        
        def getLongtimeMaxBond(self):
            return max(self.getLongtimeMaxBondLocal())
