import logging
import types, new
from espresso import esmpi as mpi

##################################################
## EXEC
##################################################
def exec_(statement) :
    """Execute arbitrary python code on all PMI workers.
    
    This allows to import modules and to define classes and functions on all workers.
    """
    __check_controller(exec_)
    # broadcast the statement
    _broadcast(_EXEC, statement)
    # locally execute the statement
    __exec_(statement)

def __exec_(statement) :
    'Internal function that executes the given statement locally.'
    # executing the statement locally
    log.info("Executing '%s'", statement)
    exec statement in globals()

# import_ is an alias for exec_
import_=exec_

##################################################
## CREATE
##################################################
def create(theClass, *args) :
    """Create an object on all workers.

    theClass describes the class of that should be instantiated.
    *args are the arguments to the constructor of the class.

    Example:
    # instantiate an object of the class with the passed name
    pmi.create("HelloWorld")
    # instantiate an object of the same class as the passed class
    pmi.create(HelloWorld)

    Note, that you can use only those classes that are know to PMI
    when this function is called, i.e. classes in modules that have
    been imported via pmi.exec_() or pmi.import_().
    """
    if isinstance(theClass, types.StringTypes) :
        theClass = eval(theClass)
    elif type(theClass) == types.TypeType :
        pass
    elif type(theClass) == types.ClassType :
        raise TypeError("pmi can't handle old-style classes. \
Please create old style classes via their names.")
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
    targs=map(__translate_arg, args)
    # broadcast creation to the workers
    _broadcast(_CREATE, theClass, oid, *targs)
    # store the destroyer object in the instance
    obj.__pmi_destroyer = __Destroyer(oid)
    return obj

def __create(theClass, oid, *args) :
    # backtranslate the arguments
    btargs=map(__backtranslate_arg, args)
    # create the new object
    obj = theClass(*btargs)
    # store the new object
    if oid in OBJECT_CACHE :
        raise InternalError("Object with oid %d is already in OBJECT_CACHE!" % oid)
    OBJECT_CACHE[oid] = obj
    log.info('Created: %s%s [oid=%d]', theClass.__name__, tuple(args), oid)


##################################################
## CLONE
##################################################
# If a class is picklable, a living instance can be cloned

##################################################
## INVOKE
##################################################
def invoke(arg0, *args) :
    """Invoke a function on all workers.

    function is the function that is to be called, *args are the
    arguments to the function. 
    Example:
    lj = _LennardJones(epsilon=1.0, sigma=1.0, cutoff=2.0)
    r = 1.2
    print(pmi.invoke(lj.computeEnergy, r))
    # equivalent:
    print(pmi.invoke(_LennardJones.computeEnergy, lj, r))

    Note, that you can use only functions that are know to PMI when
    this function is called, i.e. functions in modules that have 
    been imported via pmi.importModule().
    """
    __check_controller(invoke)
    if isinstance(arg0, types.StringTypes) :
        function = arg0
    elif isinstance(arg0, types.FunctionType) :
        function = '.'.join((arg0.__module__, arg0.__name__))
    elif isinstance(arg0, types.MethodType) :
        function = '.'.join((arg0.im_func.__module__, arg0.im_class.__name__, arg0.im_func.__name__))
        args = (arg0.im_self,) + args
    else :
        raise ValueError("pmi.invoke expects function as first argument, but got %s" % arg0)

    targs=map(__translate_arg, args)
    _broadcast(_INVOKE, function, *targs)
    log.info("Invoking: %s%s", function, tuple(args))
    function = eval(function)
    return function(*args)

def __invoke(function, *args) :
    btargs=map(__backtranslate_arg, args)
    log.info("Invoking: %s%s", function, tuple(btargs))
    function = eval(function)
    function(*btargs)

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
        log.info("Deleting OID: [%d]", self.oid)
        _broadcast(_DELETE, self.oid)
        OIDS.remove(self.oid)

def __delete(oid) :
    """Deletes the OBJECT_CACHE reference to a PMI object."""
    obj=OBJECT_CACHE[oid]
    log.info("Deleting: %s [%d]", obj, oid)
    # Delete the entry from the cache
    del OBJECT_CACHE[oid]

##################################################
## FINISH
##################################################
def finish() :
    'Ends PMI, i.e. stops all workers'
    global isFinished
    __check_controller(finish)
    if isFinished : return
    isFinished = True
    log.info('Calling all workers to stop.')
    _broadcast(_FINISH)

def __finish() :
    log.info('Finishing worker loop.')
    raise StopIteration()

##################################################
## DUMP
##################################################
def dump() :
    'Dump the object cache of PMI.'
    __check_controller(dump)
    _broadcast(_DUMP)
    log.info("OIDS=%s", str(OIDS))

def __dump() :
    log.info("OBJECT_CACHE=%s", str(OBJECT_CACHE))

##################################################
## WORKER LOOP
##################################################
def workerLoop() :
    'Starts the main worker loop. Has to be called on the workers.'
    log.info(' Entering the worker loop.')
    # On the controller, leave immediately
    if IS_CONTROLLER : 
        log.info('Leaving the worker loop.')
        return None
    
    try :
        while 1 :
            log.debug('Waiting for next PMI command.')
            __receive()
    except StopIteration :
        pass

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


def __translate_arg(obj) :
    """Internal function that translates obj into an __OID
    object if it is a PMI object instance.

    If the object is not a PMI object, returns obj untouched.
    """
    oid = id(obj)
    if oid in OIDS :
        return __OID(oid)
    else :
        return obj

def __backtranslate_arg(obj) :
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

def __check_controller(func) :
    """Check whether we are on a controller, raises a UserError if we
    are not.
    """
    if IS_WORKER : 
        raise UserError("Can't call %s on worker!" % func.__name__)

def _broadcast(*args) :
    "Internal command used to broadcast a PMI command to all workers."
    if args[0] not in _ALLCMD :
        raise ValueError('Broadcast needs a command (one of %s) as first argument.' % args)
    log.debug("Broadcasting command: %s", args)
    mpi.broadcast(mpi.world, value=args, root=CONTROLLER)

def __receive() :
    "Internal worker command that receives the broadcasts."
    message = mpi.broadcast(mpi.world, root=CONTROLLER)
    log.debug("Received command: %s", message)
    cmd = message[0]
    if cmd not in _ALLCMD :
        raise InternalError("Received bad PMI command: '%s'" % cmd)
            # determine which function to call
    cmd_func = _ALLCMD[cmd]
    args = message[1:]
    log.debug("Calling function %s%s", cmd_func.__name__, args)
    # now call the function
    cmd_func(*args)

# Command IDs
_EXEC = 'PMIEXEC'
_CREATE = 'PMICREATE'
_INVOKE = 'PMIINVOKE'
_DELETE = 'PMIDELETE'
_FINISH = 'PMIFINISH'
_DUMP = 'PMIDUMP'

# dict that associates the command IDs with their worker commands
_ALLCMD = { _EXEC : __exec_,
            _CREATE : __create,
            _INVOKE : __invoke,
            _DELETE : __delete,
            _FINISH : __finish,
            _DUMP : __dump
            }

if IS_CONTROLLER: 
    # set that stores which oids have been PMI created
    OIDS = set()
else :
    # dict that stores the objects corresponding to an oid
    OBJECT_CACHE = {}

isFinished = False

##################################################
## MODULE BODY
##################################################
if IS_CONTROLLER :
    WORKERSTR = 'Controller'
    log = logging.getLogger('espresso.pmi.controller')
else :
    WORKERSTR = 'Worker %d' % mpi.rank
    log = logging.getLogger('espresso.pmi.worker%d' % mpi.rank)


import sys, atexit
if IS_WORKER: 
    workerLoop()
    sys.exit()
else:
    atexit.register(finish)
