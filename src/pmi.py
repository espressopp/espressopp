"""The module pmi provides functions for Parallel Method Invocation
(PMI).

PMI is intended to be used in data-parallel environments, where
several threads run in parallel and can communicate via MPI.

In PMI mode, a single thread of control (a python script that runs on
the \"controller\", i.e. the MPI root task) can invoke arbitrary
functions on all other threads (the \"workers\") in parallel via
call(), invoke() and reduce(). When the function on the workers
return, the control is returned to the controller.

This model is equivalent to the \"Fork-Join execution model\" used
e.g. in OpenMP.

PMI also allows to create parallel instances of object classes via
create(), i.e. instances that have a corresponding object instance on
all workers. call(), invoke() and reduce() can be used to call
arbitrary methods of these instances.

To allow importing of python modules and to execute arbitrary code on
all workers, exec_() can be used.

Main program
------------

On the workers, the main program of a PMI script usually consists of a
single call to the function startWorkerLoop(). On the workers, this
will start an infinite loop on the workers that waits to receive the
next PMI call, while it will immediately return on the controller. On
the workers, the loop ends only, when one of the commands
finalizeWorkers() or stopWorkerLoop() is issued on the controller. A
typical PMI main program looks like this:

EXAMPLE PMI PROGRAM:

>>> # compute 2*factorial(42) in parallel
>>> import pmi
>>>
>>> # start the worker loop
>>> # on the controller, this function returns immediately
>>> pmi.startWorkerLoop()
>>>
>>> # Do the parallel computation
>>> pmi.exec_('import math')
>>> pmi.reduce('lambda a,b: a+b', 'math.factorial', 42)
>>>
>>> # exit all workers
>>> pmi.finalizeWorkers()

Instead of using pmi.finalizeWorkers() at the end of the script, you
can call registerAtExit() anywhere else, which will cause
finalizeWorkers() to be called when the python interpreter exits.

Alternatively, it is possible to use PMI in an SPMD-like fashion,
where each call to a PMI command on the controller must be accompanied
by a corresponding call on the worker. This can be either a simple
call to receive() that accepts any PMI command, or a call to the
identical PMI command. In that case, the arguments of the call to the
PMI command on the workers are ignored. In this way, it is possible to
write SPMD scripts that profit from the PMI communication patterns.

EXAMPLE SPMD PROGRAM:

>>> # compute 2*factorial(42) in parallel
>>> import pmi
>>>
>>> pmi.exec_('import math')
>>> pmi.reduce('lambda a,b: a+b', 'math.factorial', 42)

To start the worker loop, the command startWorkerLoop() can be issued
on the workers. To stop the worker loop, stopWorkerLoop() can be
issued on the controller, which will end the worker loop without
exiting the workers. 

Controller commands
-------------------

These commands can be called in the controller script. When any of
these commands is issued on a worker during the worker loop, a
UserError is raised.

* call(), invoke(), reduce() to call functions and methods in parallel
* create() to create parallel object instances
* exec_() to execute arbitrary python code in parallel
* finalizeWorkers() to stop and exit all workers
* registerAtExit() to make sure that finalizeWorkers() is called when
  python exits on the controller
* stopWorkerLoop() to interrupt the worker loop an all workers and to
  return control to the single workers

Worker commands
---------------

These commands can be called on a worker.

* startWorkerLoop() to start the worker loop
* receive() to receive a single PMI command

Useful constants
----------------

The pmi module defines the following useful constants:
* IS_CONTROLLER is True when used on the controller, False otherwise
* IS_WORKER = not IS_CONTROLLER
* ID is the numerical Id of the thread ( = mpi.world.rank)
* CONTROLLER is the numerical Id of the controller thread (the MPI root)
* WORKERSTR is a string describing the thread ('Worker #' or 'Controller')
"""
import logging
import types
from espresso import boostmpi as mpi

##################################################
## EXEC
##################################################
def exec_(statement=None) :
    """Controller command that executes arbitrary python code on all workers.
    
    This allows to import modules and to define classes and functions
    on all workers.
    Example:

    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    """
    if __checkController(exec_) :
        if statement is None :
            raise UserError('pmi.exec_ expects exactly 1 argument on controller!')
        # broadcast the statement
        _broadcast(_EXEC, statement)
        # locally execute the statement
        return __workerExec_(statement)
    else :
        return receive(_EXEC)

def __workerExec_(statement) :
    # executing the statement locally
    log.info("Executing '%s'", statement)
    exec statement in globals()

def import_(statement) :
    """import_() is an alias for exec_()."""
    exec_(statement)

##################################################
## CREATE
##################################################
def create(theClass=None, *args) :
    """Controller command that creates an object on all workers.

    "theClass" describes the (new-style) class that should be
    instantiated.
    *args are the arguments to the constructor of the class.
    Only classes that are known to PMI can be used, that is, classes
    that have been imported to pmi via exec_() or import_().

    Example:
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create("hello.HelloWorld")
    >>> print(hw)
    MPI process #0: Hello World!
    MPI process #1: Hello World!
    ...

    # Alternative.
    # Note that in this case the class has to be imported to the
    # calling module AND via PMI.
    >>> import hello
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create(hello.HelloWorld)
    >>> print(hw)
    MPI process #0: Hello World!
    MPI process #1: Hello World!
    ...
    """
    if __checkController(create) :
        if theClass is None :
          raise UserError("pmi.create expects at least 1 argument on controller")
        elif isinstance(theClass, types.StringTypes) :
            theClass = eval(theClass)
        elif type(theClass) == types.TypeType :
            pass
        elif type(theClass) == types.ClassType :
            raise TypeError("""PMI cannot create old-style classes.
            Please create old style classes via their names.
            """)
        else :
            raise ValueError("pmi.create expects class as first argument, but got %s" % theClass)

        # create the new object
        obj = theClass(*args)
        oid = id(obj)
        if oid in OIDS :
            raise _InternalError("Object with oid %d is already in OIDS!" % oid)
        OIDS.add(oid)
        log.info('Created: %s%s [oid=%d]', theClass.__name__, tuple(args), oid)
        # translate the arguments
        targs=map(__translate, args)
        # broadcast creation to the workers
        _broadcast(_CREATE, theClass, oid, *targs)
        # store the destroyer object in the instance
        obj.__pmi_destroyer = __Destroyer(oid)
        return obj
    else :
        return receive(_CREATE)

def __workerCreate(theClass, oid, *args) :
    # backtranslate the arguments
    btargs=map(__backtranslate, args)
    # create the new object
    obj = theClass(*btargs)
    # store the new object
    if oid in OBJECT_CACHE :
        raise _InternalError("Object with oid %d is already in OBJECT_CACHE!" % oid)
    OBJECT_CACHE[oid] = obj
    log.info('Created: %s%s [oid=%d]', theClass.__name__, tuple(args), oid)
    return obj

##################################################
## CLONE
##################################################
# If a class is picklable, a living instance can be cloned

##################################################
## INVOKE
##################################################
def invoke(function=None, *args) :
    """Invoke a function on all workers, gathering the return values into a list.

    function denotes the function that is to be called, *args are the
    arguments to the function.
    invoke() returns the results of the different workers as a list.
    Only functions that are known to PMI can be used, that is functions
    that have been imported to pmi via exec_() or import_().

    Example:
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    >>> messages = pmi.invoke(hw.hello())
    >>> # alternative:
    >>> messages = pmi.invoke('hello.HelloWorld.hello', hw)
    """
    if __checkController(invoke) :
        function, args = __translateInvokeArgs(function, args)
        targs=map(__translate, args)
        _broadcast(_INVOKE, function, *targs)
        log.info("Invoking: %s%s", function, tuple(args))
        function = eval(function)
        value = function(*args)
        return mpi.world.gather(value, root=CONTROLLER)
    else :
        return receive(_INVOKE)

def __workerInvoke(function, *args) :
    btargs=map(__backtranslate, args)
    log.info("Invoking: %s%s", function, tuple(btargs))
    function = eval(function)
    value = function(*btargs)
    return mpi.world.gather(value, root=CONTROLLER) 


##################################################
## CALL (INVOKE WITHOUT RESULT)
##################################################
def call(procedure=None, *args) :
    """Call a procedure on all workers, ignoring any return values.

    procedure denotes the procedure that is to be called, *args are the
    arguments to the procedure.
    If the procedure should return any results, they will be silently
    ignored.
    Only procedures that are known to PMI can be used, that is procedures
    that have been imported to pmi via exec_() or import_().
    
    Example:
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    >>> pmi.call(hw.hello)
    >>> # equivalent:
    >>> pmi.call('hello.HelloWorld', hw)
    
    Note, that you can use only procedures that are know to PMI when
    call() is called, i.e. procedures in modules that have 
    been imported via exec_().
    """
    if __checkController(call) :
        procedure, args = __translateInvokeArgs(procedure, args)
        targs=map(__translate, args)
        _broadcast(_CALL, procedure, *targs)
        log.info("Calling: %s%s", procedure, tuple(args))
        procedure = eval(procedure)
        procedure(*args)
    else :
        receive(_CALL)

def __workerCall(procedure, *args) :
    btargs=map(__backtranslate, args)
    log.info("Calling: %s%s", procedure, tuple(btargs))
    procedure = eval(procedure)
    procedure(*btargs)

def invokeNoResult(procedure, *args) :
    """invokeNoResult() is an alias for call()."""
    call(procedure, args)

##################################################
## REDUCE (INVOKE WITH REDUCED RESULT)
##################################################
def reduce(reduceOp=None, function=None, *args) :
    """Invoke a function on all workers, reducing the return values to
    a single value.

    reduceOp is the (associative) operator that is used to process the
    return values, function denotes the function that is to be called,
    *args are the arguments to the function.
    reduce() reduces the results of the different workers into a
    single value via the operation reduceOp. reduceOp is assumed to be
    associative.
    Both reduceOp and function have to be known to PMI, that is they
    must have been imported to pmi via exec_() or import_().

    Example:
    >>> pmi.exec_('import hello')
    >>> pmi.exec_('joinstr=lambda a,b: \"\\n\".join(a,b)')
    >>> hw = pmi.create('hello.HelloWorld')
    >>> print(pmi.reduce('joinstr', hw.hello()))
    >>> # equivalent:
    >>> print(
    ...   pmi.reduce('lambda a,b: \"\\n\".join(a,b)',
    ...             'hello.HelloWorld.hello', hw)
    ...             )
    """
    if __checkController(reduce) :
        # handle reduceOp argument
        if reduceOp is None :
            raise UserError("pmi.reduce expects reduceOp argument on controller")
        if isinstance(reduceOp, types.StringTypes) :
            pass
        elif isinstance(reduceOp, types.FunctionType) :
            reduceOp = '.'.join((reduceOp.__module__, reduceOp.__name__))
        else :
            raise ValueError("pmi.reduce expects function as first argument, but got %s instead" % reduceOp)
        function, args = __translateInvokeArgs(function, args)
        targs=map(__translate, args)
        _broadcast(_REDUCE, reduceOp, function, *targs)
        log.info("Reducing: %s%s", function, tuple(args))
        function = eval(function)
        reduceOp = eval(reduceOp)
        value = function(*args)
        return mpi.world.reduce(op=reduceOp, value=value,
                                root=CONTROLLER)
    else :
        return receive(_REDUCE)

def __workerReduce(reduceOp, function, *args) :
    btargs=map(__backtranslate, args)
    log.info("Invoking: %s%s", function, tuple(btargs))
    function = eval(function)
    value = function(*btargs)
    log.info("Reducing results via %s", reduceOp)
    reduceOp = eval(reduceOp)
    return mpi.world.reduce(op=reduceOp, value=value, root=CONTROLLER) 

##################################################
## AUTOMATIC OBJECT DELETION
##################################################
class __Destroyer(object) :
    """Internal class that holds the oid of a PMI instance.
    The PMI instance holds a reference to the __Destroyer object. When
    the PMI instance dies, the Destroyer also dies and removes the oid
    from OIDS."""
    def __init__(self, oid) :
        self.oid = oid
    def __del__(self) :
        if self.oid not in OIDS :
            raise _InternalError('OID %d is not in OIDS!' % self.oid)
        log.info("Adding OID to DELETED_OIDS: [%d]", self.oid)
        DELETED_OIDS.add(self.oid)
        OIDS.remove(self.oid)

def __workerDelete(*args) :
    """Deletes the OBJECT_CACHE reference to a PMI object."""
    log.info("Deleting objects: %s", args)
    for oid in args :
        obj=OBJECT_CACHE[oid]
        log.debug("  %s [%d]", obj, oid)
        # Delete the entry from the cache
        del OBJECT_CACHE[oid]

##################################################
## DUMP
##################################################
def dump() :
    """Controller function that dumps the object cache of PMI. For
    debugging purposes."""
    if __checkController(dump) :
        _broadcast(_DUMP)
        log.info("OIDS=%s", str(OIDS))
    else :
        receive(_DUMP)

def __workerDump() :
    log.info("OBJECT_CACHE=%s", str(OBJECT_CACHE))

##################################################
## WORKER LOOP CONTROL
##################################################
def registerAtExit() :
    """Controller command that registers the function
    finalizeWorkers() via atexit. 
    """
    if __checkController(registerAtExit) :
        import atexit
        atexit.register(finalizeWorkers)

def finalizeWorkers():
    """Controller command that stops and exits all workers.
    """
    stopWorkerLoop(doExit=True)

def stopWorkerLoop(doExit=False) :
    """Controller command that stops all workers.

    If \"exit\" is set, the workers exit afterwards.
    """
    if __checkController(stopWorkerLoop) :
        log.info('Calling all workers to stop.')
        _broadcast(_STOP, doExit)
    else :
        receive(_STOP)

def __workerStop(doExit) :
    import sys
    if doExit :
        log.info('Stopping worker loop and exiting worker thread.')
        sys.exit()
    else :
        log.info('Stopping worker loop.')
        raise StopIteration()

def receive(expected=None) :
    """Worker command that receives and handles the next PMI command.

    This function waits to receive and handle a single PMI command. If
    expected is not None and the received command does not equal
    \"expected\", raise a UserError.
    """
    __checkWorker(receive)
    log.debug('Waiting for next PMI command.')
    message = mpi.world.broadcast(root=CONTROLLER)
    log.debug("Received command: %s", message)
    if not hasattr(message, '__getitem__'):
        raise UserError("Received an MPI message that is not a PMI command: '%s'" % message)
    try :
        cmd = message[0]
        args = message[1:]
    except KeyError:
        raise UserError("Received an MPI message that contains an associative map: '%s'" % message)
    if cmd not in _ALLCMD :
        raise UserError("Received a message that does not start with a PMI command: '%s'" % cmd)
    elif cmd == _DELETE :
        # if delete is sent, delete the objects
        __workerDelete(*args)
        return receive(expected)
    elif expected is not None and cmd != expected :
        # otherwise test whether the command is expected
        raise UserError("Received PMI command %s but expected %s" % (cmd, expected))
    # determine which function to call
    cmd_func = _CMDS[cmd]
    log.debug("Calling function %s%s", cmd_func.__name__, args)
    # if the command is a delete, call receive once more
    if cmd == _DELETE :
        cmd_func(*args)
    else :
        return cmd_func(*args)

def startWorkerLoop() :
    """Worker command that starts the main worker loop.

    This function starts a loop that expects to receive PMI commands
    until stopWorkerLoop() or finalizeWorkers() is called on the
    controller.
    """
    # On the controller, leave immediately
    if IS_CONTROLLER :
        log.info('Entering and leaving the worker loop')
        return None

    log.info('Entering the worker loop.')
    inWorkerLoop = True

    try :
        while 1 :
            receive()
    except StopIteration :
        inWorkerLoop = False

# Metaclass
# import functools

# def method_caller(method, self, *args) :
#     print('method_caller(%s, %s)' % (self, args))
#     getattr(self, method)(self, *args)

# def getargs(*args) :
#     print('arguments=%s' % args)

# class Proxy(type):
#     def __init__(self, name, bases, dict) :
#         # init the pmi object
#         self.local = create(dict['pmiclass'])
#         for method in dict['pmilocal'] :
#             newfunc = functools.partial(method_caller, method)
#             print(newfunc.func)
#             print(newfunc.args)
# #             setattr(self, method,
# #                     functools.partial(method_caller, method))
                                     
#         return type.__init__(self, name, bases, dict)

##################################################
## CONSTANTS AND EXCEPTIONS
##################################################
# MPI task that runs the controller
CONTROLLER = 0
# whether this is a worker or a controller
IS_CONTROLLER = mpi.rank == CONTROLLER
IS_WORKER = not IS_CONTROLLER

class _InternalError(Exception) :
    """Raised when PMI has encountered an internal error.

    Hopefully, this exceptions is never raised."""
    def __str__(self) :
        return '%s: %s' % (WORKERSTR, str(self.args))
    def __repr__(self) :
        return str(self.args)

class UserError(Exception) :
    """Raised when PMI has encountered a user error.
    """
    def __str__(self) :
        return '%s: %s' % (WORKERSTR, str(self.args))
    def __repr__(self) :
        return str(self.args)

##################################################
## INTERNAL FUNTIONS
##################################################
class __OID(object) :
    """Internal class that can be sent via MPI and represents a PMI
    object instance."""
    def __init__(self, oid) :
        self.oid = oid
    def getinitargs(self) :
        return self.oid
    def __str__(self) :
        return '[%d]' % self.oid

def __checkController(func) :
    """Checks whether we are on the controller, raises a UserError if we are not.
    """
    if IS_CONTROLLER:
        return True
    else:
        if not inWorkerLoop:
            return False
        else:
            raise UserError("Cannot call %s on worker while in worker loop!" % func.__name__)

def __checkWorker(func) :
    """Checks whether we are on a worker, raises a UserError if we are not.
    """
    if IS_CONTROLLER : 
        raise UserError("Cannot call %s on the controller!" % func.__name__)

def __translate(obj) :
    """Internal function that translates obj into an __OID
    object if it is a PMI object instance.

    If the object is not a PMI object, returns obj untouched.
    """
    oid = id(obj)
    if oid in OIDS :
        return __OID(oid)
    else :
        return obj

def __backtranslate(obj) :
    """Internal worker function that backtranslates an __OID object
    into the corresponding PMI worker instance.

    If the object is not an __OID object, returns the object untouched.
    """
    if (type(obj) == __OID) :
        if obj.oid not in OBJECT_CACHE :
            raise _InternalError("OID %d was not in OBJECT_CACHE!" % obj.oid)
        return OBJECT_CACHE[obj.oid]
    else :
        return obj

def __translateInvokeArgs(arg0, args) :
    """Internal controller function that normalizes the function
    argument to invoke(), call() or reduce().
    """
    if arg0 is None :
        raise UserError("pmi.__invoke expects function argument on controller")
    if isinstance(arg0, types.StringTypes) :
        function = arg0
    elif isinstance(arg0, types.FunctionType) :
        function = '.'.join((arg0.__module__, arg0.__name__))
    elif isinstance(arg0, types.MethodType) :
        function = '.'.join((arg0.im_func.__module__, arg0.im_class.__name__, arg0.im_func.__name__))
        args = (arg0.im_self,) + args
    else :
        raise ValueError("pmi.__invoke expects function as first argument, but got %s" % arg0)
    return (function, args)

def _broadcast(*args) :
    """Internal controller function that actually broadcasts the PMI
    message.

    The function first checks whether the command is a valid PMI
    command, then it checks whether any objects have to be deleted
    before the command is broadcast, and finally it broadcasts the
    command itself.
    """
    cmd = args[0]
    if cmd not in _ALLCMD :
        raise _InternalError('_broadcast needs a command (one of %s) as first argument. Got %s instead' % (_ALLCMD, args))
    if len(DELETED_OIDS) > 0 :
        log.debug("Got %d objects in DELETED_OIDS.", len(DELETED_OIDS))
        deleteCmd = (_DELETE,) + tuple(DELETED_OIDS)
        log.debug("Broadcasting command: %s", deleteCmd)
        mpi.world.broadcast(value=deleteCmd, root=CONTROLLER)
        DELETED_OIDS.clear()
    log.debug("Broadcasting command: %s", args)
    mpi.broadcast(mpi.world, value=args, root=CONTROLLER)

# Command IDs
_ANY = 'PMIANY'
_EXEC = 'PMIEXEC'
_CREATE = 'PMICREATE'
_INVOKE = 'PMIINVOKE'
_CALL = 'PMICALL'
_REDUCE = 'PMIREDUCE'
_DELETE = 'PMIDELETE'
_STOP = 'PMISTOP'
_DUMP = 'PMIDUMP'

# dict that associates the command IDs with their worker commands
_CMDS = { _EXEC : __workerExec_,
            _CREATE : __workerCreate,
            _INVOKE : __workerInvoke,
            _CALL : __workerCall,
            _REDUCE : __workerReduce,
            _DELETE : __workerDelete,
            _STOP : __workerStop,
            _DUMP : __workerDump
            }
_ALLCMD = _CMDS.keys()

if IS_CONTROLLER: 
    # set that stores which oids have been PMI created
    OIDS = set()
    # set that stores which oids have been deleted
    DELETED_OIDS = set()
else :
    # dict that stores the objects corresponding to an oid
    OBJECT_CACHE = {}

inWorkerLoop = False
ID = mpi.rank

##################################################
## MODULE BODY
##################################################
if IS_CONTROLLER :
    WORKERSTR = 'Controller'
    log = logging.getLogger('%s.controller' % __name__)
else :
    WORKERSTR = 'Worker %d' % mpi.rank
    log = logging.getLogger('%s.worker%d' % (__name__, ID))

