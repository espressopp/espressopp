###########################################################################
#                                                                         #
#    Author: Thomas Brandes, SCAI Fraunhofer                              #
#                                                                         #
###########################################################################

from PMI import *

##########################################################################
#                                                                        #
#  Worker classes                                                        #
#                                                                        #
#    - worker class are defined directly via the Python/C++ interface    #
#                                                                        #
##########################################################################

from _Simulation import *

# nothing to be added

##########################################################################
#                                                                        #
#  Controller classes                                                    #
#                                                                        #
#    - must be inherited from class Proxy                                #
#    - all methods have the same syntax as the methods                   #
#      of the corresponding worker classes                               #
#    - each method call results in a broadcast of the arguments          #
#      to the workers that execute the command for their objects         #
#    - results of the worker methods can be gathered / reduced           #
#                                                                        #
#  Attention:                                                            #
#                                                                        #
#    - definition of the Proxy class is very error prone                 #
#    - think about doing it automatically                                #
#                                                                        #
##########################################################################

class cSimulation(Proxy):

    def __init__(self):
        Proxy.__init__(self, Simulation)
        pmiDoCommand("__init__", None, self)

    def addGroup(self, group):
        pmiDoCommand("addGroup", None, self, group)

    def getGroup(self, index):
        result = Proxy(Group)
        fresult = pmiDoCommand("getGroup", result, self, index)
        fresult.__class__ = cGroup
        return fresult

    def numberGroups(self):
        return pmiDoCommand("numberGroups", None, self)

class cGroup(Proxy):

    def __init__(self, name):
        Proxy.__init__(self, Group)
        pmiDoCommand("__init__", None, self, name)

    def getName(self):
	return pmiDoCommand("getName", None, self)

