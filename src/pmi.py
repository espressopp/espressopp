## \package PMI

##########################################################################
#                                                                        #
#   Python PMI : Parallel Method Invocation                              #
#                                                                        #
#      - This module implements parallel method invocation               #
#        where a controller process invokes methods on each worker       #
#                                                                        #
#   Comparison with PMI in C++                                           #
#                                                                        #
#     + implementation rather clear and straightforward                  #
#     + very simple to make worker classes available at controller site  #
#     + simple to combine with an SPMD model                             #
#     + simple implementation of help scripts on worker site             #
#       (parallel analysis methods)                                      #
#                                                                        #
#     - relies on Boost.Python MPI bindings                              #
#     - worker processes need Boost MPI libraries                        #
#                                                                        #
#   Open Issue: no comparison of efficiency yet                          #
#                                                                        #
#   Author: Thomas Brandes, SCAI Fraunhofer                              #
#                                                                        #
#   This file is part of the ESPResSo++ distribution                     #
#   See (http://www.espresso-pp.de)                                      #
#                                                                        #
##########################################################################

# ## This routine adds the path to the Boost MPI library that is used
# #  for the implementation of PMI.
# def setPath() :

#    import os
#    import sys

#    boostMPI = os.environ.get('BOOST_HOME')

#    if boostMPI and os.path.isdir(boostMPI):

#       sys.path.append (boostMPI + '/lib')


# ##########################################################################

# setPath()

import mpi

##########################################################################
#   pmiDoCommand (*args)                                                 #
##########################################################################

##  This routine is executed by the controller and takes care that a
#   command is executed by all workers
#
#   \param args = (methodName, redFn, proxyObject, arg1, ..., argn)
#
#   - proxyObject must be an object that has an instance on each worker
#   - methodName is the name of the method called for the instances on each worker
#   - arg1, ..., argn are the arguments for the method 
#   - redFn specifies how the result values are combined to a final result
#   - redFn == None can be used for methods without any results
#   - redFn == None can also be used if the method has the same result on each worker 
#
#   \b Examples:
#   \code
#   Proxy.__init__(obj, system)
#   pmiDoCommand("__init__", None, obj)
#   x = pmiDoCommand("get", None, obj, 0)
#   \endcode
#
#   \a This routine works also fine for constructors (methodName == "__init__"). 

def pmiDoCommand (*args):

    # in the argument list we could replace objects derived
    # from the ProxyClass by an instance of this class
    # to avoid the broadcasting of member variables

    newArgs = list (args)

    root = 0
    args  = mpi.broadcast(mpi.world, newArgs, root)

    # the controller is also a worker and will execute the command

    result = execCommand(newArgs)

    # might be that the command returns a (worker) proxy for which
    # the controller has to find its own proxy that might have more 

    return result

##########################################################################
#  reduction functions for MPI                                           #
##########################################################################

## Example reduction function that can be used for pmiDoCommand

def add(x, y):
    return x + y

##########################################################################
#                                                                        #
#   pmiWorkerLoop()                                                      #
#                                                                        #
#      - this routine is called by all workers but the controller        #
#      - it is an infinite loop waiting for new commands                 #
#      - the loop will be terminated by an empty command                 #
#                                                                        #
##########################################################################

## This routine is called by all workers waiting for new commands.
#  It is an infinite loop waiting for new commands to be executed by 
#  the worker. An empty command will terminate the loop and this routine.

def pmiWorkerLoop():

    while True :

#   receive my command from the controller

       root = 0
       command = mpi.broadcast(mpi.world, None, root)

       if (len(command) == 0):

#   empty command stands for termination of worker

          break

       else:

#   execute the command

          try:
             execCommand(command)
          except StopIteration:             
#            this is not serious as first node will pass it back to controller
             pass

##########################################################################
#                                                                        #
#   pmiFinish()                                                          #
#                                                                        #
#     - executed by the controller to finish pmiWorkerLoop()             #
#                                                                        #
##########################################################################

def pmiFinish():

      pmiDoCommand()

##########################################################################
#                                                                        #
#   execCommand (argList)                                                #
#                                                                        #
#     - this routine executes a command on a worker                      #
#     - it is sure that argList is a list and not a tuple                #
#     - important: all proxy objects are replaced by the worker objects  #
#     - argList = [methodName, redfn, proxyObj, arg1, ...., argn]        #
#                                                                        #
##########################################################################

## This routine will execute a command by a worker.
# 
#  \param argList = (methodName, redfn, proxyObj, arg1, ..., argn)
#

def execCommand(argList) :

   if (len(argList) == 0):
       return

#  argList = [methodName, redfn, proxyObject, arg1, ..., argN]

   proxyObject = argList[2]
   redFunction = argList[1]
   methodName  = argList[0]

   if isinstance (proxyObject, Proxy):
       # a proxy object is called
       workerClass   = proxyObject.workerClass
       proxyKey      = proxyObject.proxyKey
       isConstructor = (methodName == "__init__")
       
   elif isinstance (proxyObject, str):
       # import a module?
       # we just call a module function
       workerClass   = __import__(proxyObject)
       
       # ToDo: make sure that it is really a module
       proxyKey      = None
       isConstructor = False
       
   else:
      raise ArgumentError('execCommand, third arg not proxy object')

   #  now remap proxy Objects of class C to worker objects of class _C

   for i, arg in enumerate (argList):
      if isinstance (arg, Proxy) :
         if (i > 2) or (i == 2 and not isConstructor):
#            replace controller object with my worker object
             argList[i] = myObjectMap [arg.proxyKey]

#  do the same remapping also for the key argument list
#  for key, arg in keyArgs.iteritems():
#     if isinstance (arg, Proxy) :
#        replace controller object with my worker object
#        keyArgs[key] = Proxy.myObjectMap[arg.proxyKey]

   if isConstructor:

#     the constructor is called without proxyObject

      res = workerClass(*argList[3:]);

#     the worker object can be identified by the proxyKey

      myObjectMap [proxyKey] = res     

#     And for a worker object we point back to the proxy

      res.proxy = proxyObject

   elif (methodName == "delete"):

      print 'delete worker object ', argList[2]
      del(argList[2])

   else:

#     print 'proc ', mpi.rank, ' calls method ', methodName

      if proxyKey == None:
         print 'PMI: calling a module method ', methodName
         print 'module = ', workerClass
         callMethod = getattr (workerClass, methodName)
         print 'callMethod = ', callMethod
         result = callMethod(*argList[3:])
      else:
         callMethod = getattr (workerClass, methodName)
         result = callMethod(*argList[2:])

#     print 'proc ', mpi.rank, ' ready method, result = ', result

      if (result is not None) :

         if redFunction is None:

            return result

         elif isinstance (redFunction, Proxy):

#           print 'function returns proxy object, should get key', redFunction.proxyKey

#           result of the function becomes a controller object
#           identify the result worker object with the controller object

            if hasattr(result, "proxy"):

#               print 'result has already proxykey ', result.proxy.proxyKey
                return result.proxy

            else:

#               redFunction becomes the proxy of the result

#               print 'result will be registered with proxy key = ', redFunction.proxyKey
#               print 'result is of class = ', result.__class__.__name__ 

                myObjectMap [redFunction.proxyKey] = result
                result.proxy = redFunction
#               set the worker class of the proxy object
                redFunction.workerClass = result.__class__
                return redFunction

         else:

            return mpi.reduce (mpi.world, result, redFunction, 0)
   
##########################################################################
#                                                                        #
#  class Proxy                                                           #
#                                                                        #
#     - each controller class is inherited from this class Proxy         #
#     - each worker class _Example will get a controller class Example   #
#     - each object of this class has a unique key                       #
#                                                                        #
##########################################################################

## Proxy is a simple class from which all controller classes will inherit.
# 
class Proxy:

    ## class variable that contains lates key given to a proxy

    globalKey = 0

    ## dictionairy that will map a proxy key to the object for which it stands

    myObjectMap = {}

    ## Constructor of a new proxy object for a given class

    def __init__(self, workerClass):
        self.proxyKey = Proxy.globalKey
        self.workerClass = workerClass
        Proxy.globalKey += 1

    ## deletion of a proxy object deletes the object it stands for

    def delete(self):
        pmiDoCommand("delete", None, self)
        del myObjectMap[self.proxyKey]
        del(self)

##########################################################################
#                                                                        #
#  PMIExecute (controllerRoutine)                                        #
#                                                                        #
##########################################################################

##  This routine is used to execute a function in PMI mode.

def PMIExecute(controllerRoutine, *args, **kargs):

   print "I am process %d of %d." % (mpi.rank, mpi.size)

   if (mpi.rank == 0) :

#       Controller executes script

        controllerRoutine(*args, **kargs)

#       and terminates PMI execution mode

        pmiFinish()

   else :

#       Workers execute the commands they get from the host

        pmiWorkerLoop()

#  clean up the simulation by deleting the proxies

   myObjectMap = {}

   print "Process %d of %d is ready." % (mpi.rank, mpi.size)

# init the myObjectMap map
myObjectMap = {}
