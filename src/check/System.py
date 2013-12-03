"""
*************************
**espresso.check.System**
*************************

"""
def System(system, what='all'):
  if (what == 'all' or what == 'rng'):
    if system.rng == None:
      print "system has no random number generator (set with: system.rng=...)"
      return False
  if (what == 'all' or what == 'bc'):
    if system.bc == None:
      print "system has no boundary condition (set with: system.bc=...)"
      return False
  return True
