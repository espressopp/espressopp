from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class VerletListAdressLocal(_espresso.VerletListAdress):
    'The (local) verlet list AdResS'

    def __init__(self, system, cutoff, pids, exclusionlist=[]):
        'Local construction of a verlet list for AdResS'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.VerletListAdress, system, cutoff, False)
            #self.cxxclass.setAtType(self, atType)
            # check for exclusions
            if (exclusionlist != []):
                # add exclusions
                for pair in exclusionlist:
                    pid1, pid2 = pair
                    self.cxxclass.exclude(self, pid1, pid2)
            # add adress particles
            for pid in pids:
                self.cxxclass.addAdrParticle(self, pid)
            # rebuild list now
            self.cxxclass.rebuild(self)
                
            
    def totalSize(self):
        'count number of pairs in VerletList, involves global reduction'
        if pmi.workerIsActive():
            return self.cxxclass.totalSize(self)

        
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


    def addAdrParticles(self, pids, rebuild=True):
        """
        Each processor takes the broadcasted atomistic particles
        and adds it to its list.
        """
        if pmi.workerIsActive():
            for pid in pids:
                self.cxxclass.addAdrParticle(self, pid)
            if rebuild:
                # rebuild list with adress particles
                self.cxxclass.rebuild(self)

    def rebuild(self):
        if pmi.workerIsActive():
            self.cxxclass.rebuild(self)

if pmi.isController:
    class VerletListAdress(object):
        __metaclass__ = pmi.Proxy
        pmiproxydefs = dict(
            cls = 'espresso.VerletListAdressLocal',
            pmiproperty = [ 'builds' ],
            pmicall = [ 'totalSize', 'exclude', 'addAdrParticles', 'rebuild' ]
            )
