"""
* general usage
* worker/controller commands
* SPMD mode
  * use commands SPMD-like
* PMI mode:
  * use workerLoop() on workers
  * use receive() commands on workers
  * use registerAtExit() on controller
"""
import logging
import types
from espresso import boostmpi as mpi

##################################################
## EXEC
##################################################
def exec_(statement) :
    """Controller command that executes arbitrary python code on all workers.
    
    This allows to import modules and to define classes and functions
    on all workers.
    Example:

    >>> pmi.exec_("import mymodule")
    >>> myclass = pmi.create("mymodule.MyClass")
    """
    if __checkController(exec_) :
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
def create(theClass, *args) :
    """Controller command that creates an object on all workers.

    "theClass" describes the (new-style) class of that should be
    instantiated.
    *args are the arguments to the constructor of the class.
    Only classes that are known to PMI can be used, that is, classes
    that have been imported to pmi via exec_() or import_().

    Example:
    # PMI import the module hello.
    >>> pmi.exec_('import hello')
    >>> pmi.create("hello.HelloWorld")

    # Alternative.
    # Note that in this case the class has to be imported to the
    # calling module.
    >>> import hello
    >>> pmi.create(hello.HelloWorld)
    """
    if __checkController(create) :
        if isinstance(theClass, types.StringTypes) :
            theClass = eval(theClass)
        elif type(theClass) == types.TypeType :
            pass
        elif type(theClass) == types.ClassType :
            raise TypeError("""PMI can't create old-style classes.
            Please create old style classes via their names.
            """)
        else :
            raise ValueError("pmi.create expects class as first argument, but got %s" % theClass)

        # create the new object
        obj = theClass(*args)
        oid = id(obj)
        if oid in OIDS :
            raise InternalError("Object with oid %d is already in OIDS!" % oid)
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
        raise InternalError("Object with oid %d is already in OBJECT_CACHE!" % oid)
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
def invoke(function, *args) :
    """Invoke a function on all workers, gathering the return values into a list.

    function denotes the function that is to be called, *args are the
    arguments to the function.
    invoke() returns the results of the different workers as a list.
    Only functions that are known to PMI can be used, that is functions
    that have been imported to pmi via exec_() or import_().

    Example:
    pmi.exec_('import hello')
    hw = pmi.create('hello.HelloWorld')
    messages = pmi.invoke(hw.hello())
    # alternative:
    messages = pmi.invoke('hello.HelloWorld.hello', hw)
    print('\\n'.join(messages))
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
def call(procedure, *args) :
    """Call a procedure on all workers, ignoring any return values.

    procedure denotes the procedure that is to be called, *args are the
    arguments to the procedure.
    If the procedure should return any results, they will be silently
    ignored.
    Only procedures that are known to PMI can be used, that is procedures
    that have been imported to pmi via exec_() or import_().
    
    Example:
    lj = _LennardJones(epsilon=1.0, sigma=1.0, cutoff=2.0)
    r = 1.2
    print(pmi.call(lj.computeEnergy, r))
    # equivalent:
    print(pmi.call(_LennardJones.computeEnergy, lj, r))

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
def reduce(reduceOp, function, *args) :
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
    pmi.exec_('import hello')
    pmi.exec_('joinstr=lambda a,b: "\\n".join(a,b)')
    hw = pmi.create('hello.HelloWorld')
    print(pmi.reduce('joinstr', hw.hello()))
    # alternative:
    print(
      pmi.reduce(lambda a,b: '\\n'.join(a,b),
                 'hello.HelloWorld.hello', hw)
                 )
    """
    if __checkController(reduce) :
        # handle reduceOp argument
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
    reduceOp = eval(reduceOp)
    value = function(*args)
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
            raise InternalError('OID %d is not in OIDS!' % self.oid)
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
    'Dump the object cache of PMI.'
    if __checkController(dump) :
        _broadcast(_DUMP)
        log.info("OIDS=%s", str(OIDS))
    else :
        receive(_DUMP)

def __workerDump() :
    log.info("OBJECT_CACHE=%s", str(OBJECT_CACHE))

##################################################
## WORKER LOOP
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
    stopWorkerLoop(exit=True)

def stopWorkerLoop(exit=False) :
    """Controller command that stops all workers.

    If "exit" is set, the workers exit afterwards.
    """
    if __checkController(stopWorkerLoop) :
        log.info('Calling all workers to stop.')
        _broadcast(_STOP, exit)
    else :
        receive(_STOP)

def __workerStop(exit) :
    import sys
    if exit :
        log.info('Stopping worker loop and exiting worker thread.')
        sys.exit()
    else :
        log.info('Stopping worker loop.')
        raise StopIteration()

def receive(expected=None) :
    """Worker command that receives and handles the next PMI command.

    This function waits to receive and handle a single PMI command. If
    the received command is not one of "expected", raise a UserError.
    """
    __checkWorker(receive)
    log.debug('Waiting for next PMI command.')
    message = mpi.world.broadcast(root=CONTROLLER)
    log.debug("Received command: %s", message)
    cmd = message[0]
    args = message[1:]
    if cmd not in _ALLCMD :
        raise InternalError("Received an MPI message that is not a PMI command: '%s'" % cmd)
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

def workerLoop() :
    """WorkerCommand that starts the main worker loop.

    This function starts a loop that expects to receive PMI commands
    until stopWorkerLoop() or finalizeWorkers() is called on the
    controller.
    """
    log.info('Entering the worker loop.')
    # On the controller, leave immediately
    if IS_CONTROLLER : return None

    try :
        while 1 :
            receive()
    except StopIteration :
        pass

# Metaclass
import functools

def method_caller(method, self, *args) :
    print('method_caller(%s, %s)' % (self, args))
    getattr(self, method)(self, *args)

def getargs(*args) :
    print('arguments=%s' % args)

class Proxy(type):
    def __init__(self, name, bases, dict) :
        # init the pmi object
        self.local = create(dict['pmiclass'])
        for method in dict['pmilocal'] :
            newfunc = functools.partial(method_caller, method)
            print(newfunc.func)
            print(newfunc.args)
#             setattr(self, method,
#                     functools.partial(method_caller, method))
                                     
        return type.__init__(self, name, bases, dict)

##################################################
## CONSTANTS AND EXCEPTIONS
##################################################
# MPI task that runs the controller
CONTROLLER = 0
# whether this is a worker or a controller
IS_CONTROLLER = mpi.rank == CONTROLLER
IS_WORKER = not IS_CONTROLLER

class InternalError(Exception) :
    def __str__(self) :
        return '%s: %s' % (WORKERSTR, str(self.args))
    def __repr__(self) :
        return str(self.args)

class UserError(Exception) :
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

def __checkController(func) :
    """Checks whether we are on the controller, raises a UserError if we are not.
    """
    if IS_WORKER :
        if SPMD :
            return False
        else :
            raise UserError("Can't call %s on worker when not in SPMD mode!" % func.__name__)
    return True

def __checkWorker(func) :
    """Checks whether we are on a worker, raises a UserError if we are not.
    """
    if IS_CONTROLLER : 
        raise UserError("Can't call %s on controller!" % func.__name__)

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
            raise InternalError("OID %d was not in OBJECT_CACHE!" % obj.oid)
        return OBJECT_CACHE[obj.oid]
    else :
        return obj

def __translateInvokeArgs(arg0, args) :
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
    cmd = args[0]
    if cmd not in _ALLCMD :
        raise ValueError('Broadcast needs a command (one of %s) as first argument. Got %s instead' % (_ALLCMD, args))
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

# whether or not SPMD mode is used
SPMD = False

##################################################
## MODULE BODY
##################################################
if IS_CONTROLLER :
    WORKERSTR = 'Controller'
    log = logging.getLogger('%s.controller' % __name__)
else :
    WORKERSTR = 'Worker %d' % mpi.rank
    log = logging.getLogger('%s.worker%d' % (__name__, mpi.rank))
