#!/usr/bin/python

###########################################################################
#                                                                         #
#  Example script for testing PMI via Python                              #
#                                                                         #
#    Author: Thomas Brandes, SCAI Fraunhofer                              #
#                                                                         #
###########################################################################

from wrapper import *

##########################################################################
#                                                                        #
#  The script routine executed by the controller                         #
#                                                                        #
##########################################################################

def controllerRoutine():

   sim = cSimulation()

   counter = 0

   for gname in ("g1", "g2", "g3", "g4", "g5"):
       g = cGroup(gname)
       g.counter = counter 
       print 'add Group ', g.getName(), ', counter = ', counter, ' key = ', g.proxyKey
       counter = counter + 1
       sim.addGroup(g)

   print 'sim has ', sim.numberGroups(), ' Groups'

   g2 = sim.getGroup(2)
   print 'group has counter', g2.counter

   print 'Proxies: ', Proxy.myObjectMap

   del(sim)

##########################################################################
#                                                                        #
#  Main program                                                          #
#                                                                        #
##########################################################################

PMIExecute(controllerRoutine)
