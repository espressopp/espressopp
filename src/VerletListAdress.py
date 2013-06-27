"""
************************************
**VerletListAdress** - Object
************************************

The VerletListAdress is the Verlet List to be used for AdResS or H-AdResS
simulations. When creating the VerletListAdress one has to provide the system 
and specify both cutoff for the CG interaction and adrcutoff for the atomistic
interaction. Often, it is important to set the atomistic adrcutoff much bigger
than the actual interaction's cutoff would be, since also the atomistic part of
the VerletListAdress (adrPairs) is built based on the coarse-grained particle
positions. For a much larger coarse-grained cutoff it is for example possible
to also set the atomistic cutoff on the same value as the coarse-grained one.

Furthermore, the sizes of the explicit and hybrid region have to be
provided (dEx and dHy in the example below) and the center of the atomistic
region has to be set (adrCenter). In the current implementation this results in
a resolution change along the x-direction of the box. A spherical symmetry can
be obtained by only minor code changes.

Bascially the VerListAdress provides 4 lists:
---------------------------------------------

* adrZone: A list which holds all particles in the atomistic and hybrid region
* cgZone: A list which holds all particles in the coarse-grained region
* adrPairs: A list which holds all pairs which have at least one particle in the
  adrZone, i.e. in the atomistic or hybrid region
* vlPairs: A list which holds all pairs which have both particles in the cgZone,
  i.e. in the coarse-grained region

Example - creating the VerletListAdress:

>>> vl      = espresso.VerletListAdress(system, cutoff=rc, adrcut=rc, dEx=ex_size, dHy=hy_size, adrCenter=[Lx/2, Ly/2, Lz/2])

"""

from espresso import pmi
import _espresso 
import espresso
from espresso.esutil import cxxinit

class VerletListAdressLocal(_espresso.VerletListAdress):
    'The (local) verlet list AdResS'

    def __init__(self, system, cutoff, adrcut, dEx, dHy, adrCenter=[], pids=[], exclusionlist=[]):
        'Local construction of a verlet list for AdResS'
        if pmi.workerIsActive():
            cxxinit(self, _espresso.VerletListAdress, system, cutoff, adrcut, False, dEx, dHy)
            #self.cxxclass.setAtType(self, atType)
            # check for exclusions
            if (exclusionlist != []):
                # add exclusions
                for pair in exclusionlist:
                    pid1, pid2 = pair
                    self.cxxclass.exclude(self, pid1, pid2)
            # add adress particles
            if (pids != []):
                for pid in pids:
                    self.cxxclass.addAdrParticle(self, pid)
            # set adress center
            if (adrCenter != []):
                self.cxxclass.setAdrCenter(self, adrCenter[0], adrCenter[1], adrCenter[2])
            
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
