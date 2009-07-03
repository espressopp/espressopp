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
* exec_() and import_() to execute arbitrary python code in parallel
  and to import classes and functions into the global namespace of
  pmi.
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
* call(), invoke(), reduce(), create() and exec_() to receive a single
  corresponding PMI command. Note that these commands will ignore any
  arguments when called on a worker.

PMI Proxy metaclass
-------------------

The Proxy metaclass can be used to easily generate front-end classes
to distributed PMI classes.
.
.
.

Useful constants and variables
------------------------------

The pmi module defines the following useful constants:
* IS_CONTROLLER is True when used on the controller, False otherwise
* IS_WORKER = not IS_CONTROLLER
* ID is the numerical Id of the thread ( = mpi.world.rank)
* CONTROLLER is the numerical Id of the controller thread (the MPI root)
* WORKERSTR is a string describing the thread ('Worker #' or 'Controller')

and the following variables:
* inWorkerLoop is True, if PMI currently executes the worker loop on
  the workers.
"""
import logging, types, sys, inspect
from espresso import boostmpi as mpi

##################################################
## EXEC
##################################################
def exec_(*args) :
    """Controller command that executes arbitrary python code on all workers.
    
    exec_() allows to execute arbitrary Python code on all workers.
    It can be used to import modules, or to define classes and
    functions on all workers.

    if arg is a string, the string is the Python code to be
    executed. If arg is any other object, exec_ will try to get its
    source code via inspect.getsource() and execute it.
    Note, that this means that if you apply exec_() to a function,
    that the function will be defined on the workers, but it will not
    execute the function itself. Note also, that function or class
    definitions will affect the top-level pmi namespace and not be
    part of the module they are defined in.

    Example:
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')

    Example2:
    >>> class HelloWorld(object):
    >>>   def hello():
    >>>     print('Hello World!')
    >>>
    >>> pmi.exec_(HelloWorld)
    >>> hw = pmi.create(HelloWorld)
    """
    if __checkController(exec_) :
        if len(args) == 0:
            raise UserError('pmi.exec_ expects exactly 1 argument on controller!')
            
        # broadcast the statement
        _broadcast(_EXEC, *args)
        # locally execute the statement
        return __workerExec_(*args)
    else :
        return receive(_EXEC)

def __workerExec_(*statements) :
    # executing the statement locally
    for statement in statements:
        log.info("Executing '%s'", statement)
        exec statement in globals()

def import_(statement) :
    """import_() is an alias for exec_()."""
    exec_(statement)

##################################################
## EXECFILE
##################################################
def execfile_(file):
    if __checkController(execfile_):
        _broadcast(_EXECFILE, file)
        return __workerExecfile_(file)
    else:
        return receive(_EXECFILE)

def __workerExecfile_(file):
    log.info("Executing file '%s'", file)
    execfile(file, globals())

##################################################
## CREATE
##################################################
def create(cls=None, *args, **kwds) :
    """Controller command that creates an object on all workers.

    cls describes the (new-style) class that should be
    instantiated.
    *args are the arguments to the constructor of the class.
    Only classes that are known to PMI can be used, that is, classes
    that have been imported to pmi via exec_() or import_().

    Example:
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
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
        cls = _translateClass(cls)

        # generate a new oid
        oid = __OID()

        # translate the arguments
        ckwds, targs, tkwds = __translateArgs(args, kwds)
        # broadcast creation to the workers
        _broadcast(_CREATE, cls, oid, *targs, **tkwds)

        log.info('Creating: %s [%s]', 
                 __formatCall(cls.__name__, args, ckwds), oid)

        # create the instance
        obj = cls(*args, **ckwds)
        # store the oid in the instance
        obj.__pmioid = oid
        obj.__pmidestroyer = __Destroyer(oid)

        return obj
    else :
        return receive(_CREATE)

def __workerCreate(cls, oid, *targs, **tkwds) :
    # backtranslate the arguments
    args, kwds = __backtranslateArgs(targs, tkwds)
    log.info('Creating: %s [%s]'
             % (__formatCall(cls.__name__, args, kwds), oid))
    # create the new object
    obj = cls(*args, **kwds)
    # store the new object
    if oid in OBJECT_CACHE :
        raise _InternalError("Object [%s] is already in OBJECT_CACHE!" % oid)
    OBJECT_CACHE[oid] = obj
    return obj

##################################################
## CLONE
##################################################
# If a class is picklable, a living instance can be cloned

##################################################
## INVOKE
##################################################
def invoke(function=None, *args, **kwds) :
    """Invoke a function on all workers, gathering the return values into a list.

    function denotes the function that is to be called, *args and
    **kwds are the arguments to the function. If **kwds contains keys
    that start with with the prefix '__pmictr_', they are stripped of
    the prefix and are passed only to the controller.

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
        function = __translateInvokeArg(function)
        ckwds, targs, tkwds = __translateArgs(args, kwds)
        _broadcast(_INVOKE, function, *targs, **tkwds)

        log.info("Invoking: %s", __formatCall(function, args, ckwds))
        value = __invoke(function, args, ckwds)
        return mpi.world.gather(value, root=CONTROLLER)
    else :
        return receive(_INVOKE)

def __workerInvoke(function, *targs, **tkwds) :
    args, kwds = __backtranslateArgs(targs, tkwds)
    log.info("Invoking: %s", __formatCall(function, args, kwds))
    value = __invoke(function, args, kwds)
    return mpi.world.gather(value, root=CONTROLLER) 

##################################################
## CALL (INVOKE WITHOUT RESULT)
##################################################
def call(function=None, *args, **kwds) :
    """Call a function on all workers, returning only the return value on the controller.

    function denotes the function that is to be called, *args and
    **kwds are the arguments to the function. If **kwds contains keys
    that start with with the prefix '__pmictr_', they are stripped of
    the prefix and are passed only to the controller.
    If the function should return any results on the workers, they
    will be silently ignored. Only the result on the controller will
    be returned by call().
    Only functions that are known to PMI can be used, that is functions
    that have been imported to pmi via exec_() or import_().
    
    Example:
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    >>> pmi.call(hw.hello)
    >>> # equivalent:
    >>> pmi.call('hello.HelloWorld', hw)
    
    Note, that you can use only functions that are know to PMI when
    call() is called, i.e. functions in modules that have 
    been imported via exec_().
    """
    if __checkController(call) :
        function = __translateInvokeArg(function)
        ckwds, targs, tkwds = __translateArgs(args, kwds)
        _broadcast(_CALL, function, *targs, **tkwds)
        log.info("Calling: %s", __formatCall(function, args, ckwds))
        return __invoke(function, args, ckwds)
    else :
        receive(_CALL)

def __workerCall(function, *targs, **tkwds) :
    args, kwds = __backtranslateArgs(targs, tkwds)
    log.info("Calling: %s", __formatCall(function, args, kwds))
    __invoke(function, args, kwds)

def invokeNoResult(function, *args, **kwds) :
    """invokeNoResult() is an alias for call()."""
    call(function, *args, **kwds)

##################################################
## REDUCE (INVOKE WITH REDUCED RESULT)
##################################################
def reduce(reduceOp=None, function=None, *args, **kwds) :
    """Invoke a function on all workers, reducing the return values to
    a single value.

    reduceOp is the (associative) operator that is used to process the
    return values, function denotes the function that is to be called,
    *args and **kwds are the arguments to the function. If **kwds
    contains keys that start with with the prefix '__pmictr_', they
    are stripped of the prefix and are passed only to the controller.

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
        function = __translateInvokeArg(function)
        ckwds, targs, tkwds = __translateArgs(args, kwds)
        _broadcast(_REDUCE, reduceOp, function, *targs, **tkwds)
        log.info("Reducing: %s", __formatCall(function, args, ckwds))
        value = __invoke(function, args, ckwds)
        log.info("Reducing results via %s", reduceOp)
        reduceOp = eval(reduceOp)
        return mpi.world.reduce(op=reduceOp, value=value, root=CONTROLLER)
    else :
        return receive(_REDUCE)

def __workerReduce(reduceOp, function, *targs, **tkwds) :
    args, kwds = __backtranslateArgs(targs, tkwds)
    log.info("Reducing: %s", __formatCall(function, args, kwds))
    value = __invoke(function, args, kwds)
    log.info("Reducing results via %s", reduceOp)
    reduceOp = eval(reduceOp)
    return mpi.world.reduce(op=reduceOp, value=value, root=CONTROLLER) 

##################################################
## AUTOMATIC OBJECT DELETION
##################################################
def __delete() :
    """Internal controller command that deletes PMI objects on the
    workers that have already been deleted on the controller.
    """
    global DELETED_OIDS
    if len(DELETED_OIDS) > 0 :
        log.debug("Got %d objects in DELETED_OIDS.", len(DELETED_OIDS))
        __broadcastCmd(_DELETE, *DELETED_OIDS)
        DELETED_OIDS = []

def __workerDelete(*args) :
    """Deletes the OBJECT_CACHE reference to a PMI object."""
    log.info("Deleting oids: %s", args)
    for oid in args :
        obj=OBJECT_CACHE[oid]
        log.debug("  %s [%s]" % (obj, oid))
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
    else :
        receive(_DUMP)

def __workerDump() :
    import pprint
    log.info("OBJECT_CACHE=%s", pprint.pformat(OBJECT_CACHE))

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
    else:
        raise UserError('Cannot call registerAtExit on worker!')

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
        raise UserError('Cannot call stopWorkerLoop on worker!')

def __workerStop(doExit) :
    if doExit :
        log.info('Stopping worker loop and exiting worker thread.')
        sys.exit()
    else :
        log.info('Stopping worker loop.')
        raise StopIteration()

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

##################################################
## PROXY METACLASS
##################################################

class Proxy(type):
    class Initializer(object):
        """Creates the PMI subject of a class.
        """
        def __init__(self, cls, oldinit):
            self.cls = cls
            if oldinit is not None:
                self.oldinit = oldinit
        def __call__(self, method_self, *args, **kwds):
            # create the pmi subject
            log.info('PMI.Proxy class %s is creating pmi subject of type %s',
                     self.cls, self.cls.pmiobjectclassdef)
            self.cls._normalizeClass()
            method_self.pmiobject = create(self.cls.pmiobjectclass, *args, **kwds)
            if hasattr(self, 'oldinit'):
                self.oldinit(method_self)
#                self.oldinit(method_self, *args, **kwds)

    class LocalCaller(object):
        def __init__(self, cls, methodName):
            self.cls = cls
            self.methodName = methodName
        def __call__(self, method_self, *args, **kwds):
            self.cls._normalizeClass()
            method = getattr(self.cls.pmiobjectclass, self.methodName)
            return method(method_self.pmiobject, *args, **kwds)

    class PMICaller(object):
        def __init__(self, cls, methodName):
            self.cls = cls
            self.methodName = methodName
        def __call__(self, method_self, *args, **kwds):
            self.cls._normalizeClass()
            method = getattr(self.cls.pmiobjectclass, self.methodName)
            return call(method, method_self.pmiobject, *args, **kwds)

    class PMIInvoker(object):
        def __init__(self, cls, methodName):
            self.cls = cls
            self.methodName = methodName
        def __call__(self, method_self, *args, **kwds):
            self.cls._normalizeClass()
            method = getattr(self.cls.pmiobjectclass, self.methodName)
            return invoke(method, method_self.pmiobject, *args, **kwds)

    class PropertyLocalGetter(object):
        def __init__(self, cls, propName):
            self.cls = cls
            self.propName = propName
        def __call__(self, method_self):
            self.cls._normalizeClass()
            property = getattr(self.cls.pmiobjectclass, self.propName)
            self.getter = getattr(property, 'fget')
            return self.getter(method_self.pmiobject)

    class PropertyPMISetter(object):
        def __init__(self, cls, propName):
            self.cls = cls
            self.propName = propName
        def __call__(self, method_self, val):
            setter = '.'.join(
                (self.cls.pmiobjectclassdef,
                 self.propName, 
                 'fset'))
            return call(setter, method_self.pmiobject, val)

    def _normalizeClass(cls):
        if not hasattr(cls, 'pmiobjectclass'):
            cls.pmiobjectclass = _translateClass(cls.pmiobjectclassdef)

    def __addMethod(cls, methodName, caller):
        newMethod = types.MethodType(caller, None, cls)
        setattr(cls, methodName, newMethod)

    def __init__(cls, name, bases, dict):
        if 'pmiproxydefs' not in dict:
            raise(UserError('PMI proxy class %s has to define pmiproxydefs.' % name))
        log.info('Making %s a proxy class.' % name)
        defs = dict['pmiproxydefs']

        if 'subjectclass' in defs:
            cls.pmiobjectclassdef = defs['subjectclass']
            if hasattr(cls, '__init__'): 
                oldinit = getattr(cls, '__init__')
            else: 
                oldinit = None
            cls.__addMethod('__init__', 
                            Proxy.Initializer(cls, oldinit))
                        
        if 'localcall' in defs:
            for methodName in defs['localcall']:
                cls.__addMethod(methodName, 
                                Proxy.LocalCaller(cls, methodName))

        if 'pmicall' in defs:
            for methodName in defs['pmicall']:
                cls.__addMethod(methodName, 
                                Proxy.PMICaller(cls, methodName))

        if 'pmiinvoke' in defs:
            for methodName in defs['pmiinvoke']:
                cls.__addMethod(methodName, 
                                Proxy.PMIInvoker(cls, methodName))

        if 'pmiproperty' in defs:
            for propname in defs['pmiproperty']:
                newprop = property(
                    Proxy.PropertyLocalGetter(cls, propname),
                    Proxy.PropertyPMISetter(cls, propname))
                setattr(cls, propname, newprop)

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
## BROADCAST AND RECEIVE
##################################################
def _broadcast(cmd, *args, **kwds) :
    """Internal controller command that actually broadcasts a PMI command.

    The function first checks whether cmd is a valid PMI command, then
    it checks whether any objects have to be deleted before the
    command is broadcast, and finally it broadcasts the command
    itself.
    """
    __delete()
    __broadcastCmd(cmd, *args, **kwds)

def __broadcastCmd(cmd, *args, **kwds) :
    if not _checkCommand(cmd) :
        raise _InternalError('_broadcast needs a PMI command as first argument. Got %s instead!' % cmd)
    cmd = __CMD(cmd, args, kwds)
    log.debug("Broadcasting command: %s", cmd)
    mpi.broadcast(mpi.world, value=cmd, root=CONTROLLER)

def receive(expected=None) :
    """Worker command that receives and handles the next PMI command.

    This function waits to receive and handle a single PMI command. If
    expected is not None and the received command does not equal
    \"expected\", raise a UserError.
    """
    __checkWorker(receive)
    if expected is None:
        log.debug('Waiting for next PMI command.')
    else:
        log.debug('Waiting for PMI command %s.', _CMD[expected][0])
    message = mpi.world.broadcast(root=CONTROLLER)
    log.debug("Received message: %s", message)
    if type(message) != __CMD:
        raise UserError("Received an MPI message that is not a PMI command: '%s'" % str(message))
    cmd = message.cmd
    args = message.args
    kwds = message.kwds
    if cmd == _DELETE :
        # if delete is sent, delete the objects
        __workerDelete(*args)
        # recursively call receive once more
        return receive(expected)
    elif expected is not None and cmd != expected :
        # otherwise test whether the command is expected
        raise UserError("Received PMI command %s but expected %s" % (_CMD[cmd][0], _CMD[expected][0]))
    # determine which function to call
    cmd_func = _CMD[cmd][1]
    log.debug("Calling function %s", __formatCall(cmd_func.__name__, args, kwds))
    return cmd_func(*args, **kwds)

##################################################
## INTERNAL FUNTIONS
##################################################
class __OID(object) :
    """Internal class that represents a PMI object id.
    
    An instance of this class can be pickled, so that it can be sent
    via MPI, and it is hashable, so that it can be used as a hash key
    (for OBJECT_CACHE).
    """
    def __init__(self) :
        self.id = id(self)
        return object.__init__(self)
    def __str__(self):
        return 'oid=0x%x' % self.id
    def __hash__(self):
        return self.id
    def __eq__(self, obj):
        return self.id == obj.id
    def __getstate__(self):
        return self.id
    def __setstate__(self, id):
        self.id = id

if IS_CONTROLLER:
    class __Destroyer(object):
        def __init__(self, oid):
            self.oid = oid
            return object.__init__(self)
        def __del__(self):
            log.info("Adding OID to DELETED_OIDS: [%s]", self.oid)
            DELETED_OIDS.append(self.oid)

class __CMD(object) :
    """Internal class that can be sent via MPI and represents a PMI command.
    """
    def __init__(self, cmd, args=None, kwds=None) :
        if not _checkCommand(cmd):
            raise _InternalError('Created __CMD object with invalid PMI command %s' % cmd)
        self.cmd = cmd
        self.args = args
        self.kwds = kwds
    def __str__(self):
        sargs = [_CMD[self.cmd][0]]
        if hasattr(self, 'args'):
            sargs.append(str(map(str, self.args)))
        if hasattr(self, 'kwds'):
            sargs.append(str(self.kwds))
        return 'PMICMD(%s)' % (', '.join(sargs))
    def __getstate__(self):
        state = (self.cmd, self.args, self.kwds)
        return state
    def __setstate__(self, state):
        self.cmd, self.args, self.kwds = state

def _checkCommand(cmd):
    return 0 <= cmd < _MAXCMD

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

def _translateClass(cls):
    """Returns the class object of the class decscribed by cls.
    """
    if cls is None :
        raise UserError("pmi.create expects at least 1 argument on controller")
    elif isinstance(cls, types.StringTypes) :
        return eval(cls)
    elif isinstance(cls, types.TypeType) :
        return cls
    elif isinstance(cls, types.ClassType) :
        raise TypeError("""PMI cannot use old-style classes.
        Please create old style classes via their names.
        """)
    else :
        raise ValueError("__translateClass expects class as argument, but got %s" % cls)


def __mapArgs(func, args, kwds):
    """Internal function that maps a function to the list args and to
    the values of the dict kwds. Used by __translateArsg and
    __backtranslateArgs.
    """
    targs = map(func, args)
    tkwds = {}
    for k, v in kwds.iteritems():
        tkwds[k] = func(v)
    return targs, tkwds
    
def _translateOID(obj) :
    """Internal function that translates obj into an __OID
    object if it is a PMI object instance.
        
    If the object is not a PMI object, returns obj untouched.
    """
    if hasattr(obj, 'pmiobject'):
        # if it is a proxy, send its pmiobjects oid
        return obj.pmiobject.__pmioid
    if hasattr(obj, '__pmioid'):
        # if it is a pmi object, send the oid
        return obj.__pmioid
    else:
        return obj

def __translateArgs(args, kwds):
    """Internal function that translates all PMI object instances that
    occur in args or kwds into __OID objects that can be sent to the
    workers.
    """
    
    workerKwds={}
    controllerKwds={}
    for k in kwds.keys():
        if k.startswith('__pmictr_'):
            knew = k[9:]
            controllerKwds[knew] = kwds[k]
        else:
            v = kwds[k]
            workerKwds[k] = v
            if k not in controllerKwds:
                controllerKwds[k] = v

    targs, tWorkerKwds = __mapArgs(_translateOID, args, workerKwds)

    return controllerKwds, targs, tWorkerKwds

def _backtranslateOID(obj) :
    """Internal worker function that backtranslates an __OID object
    into the corresponding PMI worker instance.
    
    If the object is not an __OID object, returns the object untouched.
    """
    if type(obj) is __OID:
        if obj in OBJECT_CACHE.keys():
            return OBJECT_CACHE[obj]
        else:
            raise _InternalError("Object [%s] is not in OBJECT_CACHE" % obj)
    else :
        return obj

def __backtranslateArgs(targs, tkwds):
    """Internal function that backtranslates all __OID object
    instances that occur in args or kwds into the cached PMI objects
    on the worker.
    """
    return __mapArgs(_backtranslateOID, targs, tkwds)

# Wrapper that allows to pickle a method
class __Method(object) :
    def __init__(self, funcname, im_self, cls) :
        self.__name__ = funcname
        self.im_self = im_self
        self.im_class = cls
        self.__determineMethod()
    def __getstate__(self):
        return (self.__name__, _translateOID(self.im_self), self.im_class)
    def __setstate__(self, state):
        self.__name__, self.im_self, self.im_class = state
        self.im_self = _backtranslateOID(self.im_self)
        self.__determineMethod()
    def __determineMethod(self):
        for cls in self.im_class.mro():
            if hasattr(cls, self.__name__):
                function = getattr(cls, self.__name__)
                self.method = function.__get__(self.im_self, cls)
                break
    def __call__(self, *args, **kwds):
        return self.method(*args, **kwds)

# translate arguments to invoke
def __translateInvokeArg(arg0) :
    """Internal controller function that normalizes the function
    argument to invoke(), call() or reduce().
    """
    if arg0 is None :
        raise UserError("pmi.__invoke expects function argument on controller")
    if isinstance(arg0, types.StringTypes):
        function = arg0
    elif isinstance(arg0, (types.FunctionType, types.BuiltinFunctionType)):
        if arg0.__name__ == '<lambda>':
            raise ValueError("pmi.__translateInvokeArgs cannot handle lambda functions")
        function = arg0
    elif isinstance(arg0, (types.MethodType, types.BuiltinMethodType)):
        function = __Method(arg0.im_func.__name__, arg0.im_self, arg0.im_class)
    else :
        raise ValueError("pmi.__translateInvokeArgs expects function as first argument, but got %s" % arg0)
    return function

def __invoke(function, args, kwds):
    if isinstance(function, types.StringTypes):
        function = eval(function, globals())
    log.debug("Invoking function: %s",
              __formatCall(function.__name__, args, kwds))
    return function(*args, **kwds)

def __formatCall(function, args, kwds) :
    def formatArgs(args, kwds) :
        arglist = [repr(arg) for arg in args]
        for k, v in kwds.iteritems():
            arglist.append('%s=%r' % (k, repr(v)))
        return ', '.join(arglist)

    return '%s(%s)' % (function, formatArgs(args, kwds))

# map of command names and associated worker functions
_CMD = [ ('EXEC', __workerExec_),
         ('EXECFILE', __workerExecfile_),
          ('CREATE', __workerCreate),
          ('INVOKE', __workerInvoke),
          ('CALL', __workerCall),
          ('REDUCE', __workerReduce),
          ('DELETE', __workerDelete),
          ('STOP', __workerStop),
          ('DUMP', __workerDump) ]
_MAXCMD = len(_CMD)

# define the numerical constants to be used
for i in range(len(_CMD)) :
    exec '_%s=%s' % (_CMD[i][0],i) in globals()
del i

if IS_CONTROLLER: 
    # set that stores which oids have been deleted
    DELETED_OIDS = []
else :
    # dict that stores the objects corresponding to an oid
    OBJECT_CACHE = {}

inWorkerLoop = False
ID = mpi.rank
#Globals = {}

##################################################
## MODULE BODY
##################################################
if IS_CONTROLLER :
    WORKERSTR = 'Controller'
    log = logging.getLogger('%s.controller' % __name__)
else :
    WORKERSTR = 'Worker %d' % mpi.rank
    log = logging.getLogger('%s.worker%d' % (__name__, ID))

