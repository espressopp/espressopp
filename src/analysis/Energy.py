"""
****************************
**espresso.analysis.Energy**
****************************

"""
import espresso

class EnergyPot():
    'Potential energy of the system.'
    def __init__(self, system, per_atom=False):
        self.system   = system
        self.per_atom = per_atom
        
    def compute(self):
        EPot = 0.0
        for k in range(self.system.getNumberOfInteractions()):
          EPot += self.system.getInteraction(k).computeEnergy()
        if self.per_atom:
          NPart  = espresso.analysis.NPart(self.system).compute()
          return EPot / NPart
        else:
          return EPot

class EnergyKin():
    'Kinetic energy of the system.'
    def __init__(self, system, per_atom=False):
        self.system   = system
        self.per_atom = per_atom
        
    def compute(self):
      NPart  = espresso.analysis.NPart(self.system).compute()
      T      = espresso.analysis.Temperature(self.system).compute()
      EKin   = (3.0/2.0) * NPart * T
      if self.per_atom:
          return EKin / NPart
      else:
          return EKin

class EnergyTot():
    'Total energy (EKin + EPot) of the system.'
    def __init__(self, system, per_atom=False):
        self.system   = system
        self.per_atom = per_atom
        
    def compute(self):
      NPart  = espresso.analysis.NPart(self.system).compute()
      T      = espresso.analysis.Temperature(self.system).compute()
      EKin   = (3.0/2.0) * NPart * T
      EPot   = 0.0
      for k in range(self.system.getNumberOfInteractions()):
        EPot += self.system.getInteraction(k).computeEnergy()
      if self.per_atom:
        return (EPot + EKin) / NPart
      else:
        return (EPot + EKin)



