#!/usr/bin/python

###########################################################################
#                                                                         #
#  Example script for testing the LAMMPS Python/C++ wrapper               #
#                                                                         #
#    Author: Thomas Brandes, SCAI Fraunhofer                              #
#                                                                         #
###########################################################################

from _Simulation import *

sim = Simulation()

counter = 0

for gname in ("g1", "g2", "g3", "g4", "g5"):
    g = Group(gname)
    g.counter = counter 
    counter = counter + 1
    sim.addGroup(g)

for i in range(sim.numberGroups()):
    g = sim.getGroup(i)
    if hasattr(g, "counter"):
       print 'group ', g.getName(), ', python counter = ', g.counter
    else:
       print 'group ', g.getName(), ', no python counter'

g2 = sim.getGroup(2)
print 'group has name ', g2.getName()
print 'group has counter', g2.counter

del(sim)
