from espressopp import pmi
import _espressopp 
import espressopp
from espressopp.esutil import cxxinit

class VirtualVerletListLocal(_espressopp.VirtualVerletList):
    'The (local) verlet list.'

    def __init__(self, system, cutoff, ftpl, exclusionlist=[]):
        'Local construction of a verlet list'
        if pmi.workerIsActive():
            if (exclusionlist == []):
                # rebuild list in constructor
                cxxinit(self, _espressopp.VirtualVerletList, system, cutoff, ftpl, True )
            else:
                # do not rebuild list in constructor
                cxxinit(self, _espressopp.VirtualVerletList, system, cutoff, ftpl, False)
                # add exclusions
                for pair in exclusionlist:
                    pid1, pid2 = pair
                    self.cxxclass.exclude(self, pid1, pid2)
                # now rebuild list with exclusions
                self.cxxclass.rebuild(self)
                
            
    def totalSize(self):
        'count number of pairs in VirtualVerletList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.totalSize(self)

    def setCellList(self, cl):
        'set special cellList'
        if pmi.workerIsActive():
            return self.cxxclass.setCellList(self, cl)
        
    def localSize(self):
        'count number of pairs in local VirtualVerletList'
        if pmi.workerIsActive():
            return self.cxxclass.localSize(self)
        
    def exclude(self, exclusionlist):
        """
        Each processor takes the broadcasted exclusion list
        and adds it to its list.
        """
        if pmi.workerIsActive():
            for pair in exclusionlist:
                pid1, pid2 = pair
                self.cxxclass.exclude(self, pid1, pid2)
            # rebuild list with exclusions
            self.cxxclass.rebuild(self)
            
    def getAllPairs(self):
        'return the pairs of the local verlet list'
        if pmi.workerIsActive():
            pairs=[]
            npairs=self.localSize()
            for i in range(npairs):
              pair=self.cxxclass.getPair(self, i+1)
              pairs.append(pair)
            return pairs 

    def mapTypeToVerletList(self, types, vl):
		#add a tuple of type. Atoms of that type are put in a seperate vl
		if pmi.workerIsActive():
			for t in types:
				for u in types:
					self.cxxclass.addTypeToVLMap(self, t, u, vl)

if pmi.isController:
  class VirtualVerletList(object):
    __metaclass__ = pmi.Proxy
    pmiproxydefs = dict(
      cls = 'espressopp.VirtualVerletListLocal',
      pmiproperty = [ 'builds' ],
      pmicall = [ 'totalSize', 'exclude', 'connect', 'disconnect', 'getVerletCutoff', 'setCellList','mapTypeToVerletList'],
      pmiinvoke = [ 'getAllPairs' ]
    )
