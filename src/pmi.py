#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.


# PMI - Parallel Method Invocation
# Copyright (C) 2009,2010 Olaf Lenz
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
"""
**************
espressopp.pmi
**************

Parallel Method Invocation (PMI) allows users to write serial Python scripts
that use functions and classes that are executed in parallel.

PMI is intended to be used in data-parallel environments, where
several threads run in parallel and can communicate via MPI.

In PMI mode, a single thread of control (a python script that runs on
the *controller*, i.e. the MPI root task) can invoke arbitrary
functions on all other threads (the *workers*) in parallel via
`call()`, `invoke()` and `reduce()`. When the function on the workers
return, the control is returned to the controller.

This model is equivalent to the \"Fork-Join execution model\" used
e.g. in OpenMP.

PMI also allows to create parallel instances of object classes via
`create()`, i.e. instances that have a corresponding object instance
on all workers. `call()`, `invoke()` and `reduce()` can be used to
call arbitrary methods of these instances.

to execute arbitrary code on all workers, `exec_()` can be used, and
to import python modules to all workers, use 'import_()'.

**Main program**

On the workers, the main program of a PMI script usually consists of a
single call to the function `startWorkerLoop()`. On the workers, this
will start an infinite loop on the workers that waits to receive the
next PMI call, while it will immediately return on the controller. On
the workers, the loop ends only, when one of the commands
`finalizeWorkers()` or `stopWorkerLoop()` is issued on the
controller. A typical PMI main program looks like this:

>>> # compute 2*factorial(42) in parallel
>>> import pmi
>>>
>>> # start the worker loop
>>> # on the controller, this function returns immediately
>>> pmi.startWorkerLoop()
>>>
>>> # Do the parallel computation
>>> pmi.import_('math')
>>> pmi.reduce('lambda a,b: a+b', 'math.factorial', 42)
>>>
>>> # exit all workers
>>> pmi.finalizeWorkers()

Instead of using `finalizeWorkers()` at the end of the script, you can
call `registerAtExit()` anywhere else, which will cause
`finalizeWorkers()` to be called when the python interpreter exits.

Alternatively, it is possible to use PMI in an SPMD-like fashion,
where each call to a PMI command on the controller must be accompanied
by a corresponding call on the worker. This can be either a simple
call to `receive()` that accepts any PMI command, or a call to the
identical PMI command. In that case, the arguments of the call to the
PMI command on the workers are ignored. In this way, it is possible to
write SPMD scripts that profit from the PMI communication patterns.

>>> # compute 2*factorial(42) in parallel
>>> import pmi
>>>
>>> pmi.exec_('import math')
>>> pmi.reduce('lambda a,b: a+b', 'math.factorial', 42)

To start the worker loop, the command `startWorkerLoop()` can be
issued on the workers. To stop the worker loop, `stopWorkerLoop()` can
be issued on the controller, which will end the worker loop without
exiting the workers.

**Controller commands**

These commands can be called in the controller script. When any of
these commands is issued on a worker during the worker loop, a
`UserError` is raised.

* `call()`, `invoke()`, `reduce()` to call functions and methods in parallel
* `create()` to create parallel object instances
* `exec_()` and `import_()` to execute arbitrary python code in
  parallel and to import classes and functions into the global
  namespace of pmi.
* `sync()` to make sure that all deleted PMI objects have been deleted.
* `finalizeWorkers()` to stop and exit all workers
* `registerAtExit()` to make sure that finalizeWorkers() is called when
  python exits on the controller
* `stopWorkerLoop()` to interrupt the worker loop an all workers and to
  return control to the single workers

**Worker commands**

These commands can be called on a worker.

* `startWorkerLoop()` to start the worker loop
* `receive()` to receive a single PMI command
* `call()`, `invoke()`, `reduce()`, `create()` and `exec_()` to
  receive a single corresponding PMI command. Note that these commands
  will ignore any arguments when called on a worker.

**PMI Proxy metaclass**

The `Proxy` metaclass can be used to easily generate front-end classes
to distributed PMI classes.
.
.
.

**Useful constants and variables**

The pmi module defines the following useful constants and variables:

* `isController` is True when used on the controller, False otherwise
* `isWorker` = not isController
* `ID` is the rank of the MPI task
* `CONTROLLER` is the rank of the Controller (normally the MPI root)
* `workerStr` is a string describing the thread ('Worker #' or 'Controller')
* `inWorkerLoop` is True, if PMI currently executes the worker loop on
  the workers.

"""
import logging, types, sys, inspect, os

from collections.abc import Iterable

__author__ = 'Olaf Lenz'
__email__ = 'olaf at lenz dot name'
__version__ = '1.0'
__all__ = [
    'exec_', 'import_', 'execfile_',
    'create', 'call', 'invoke', 'reduce', 'localcall',
    'sync', 'receive',
    'startWorkerLoop',
    'finalizeWorkers', 'stopWorkerLoop', 'registerAtExit',
    'Proxy',
    'rank', 'size', 'CONTROLLER',
    'isController', 'isWorker',
    'workerStr', 'inWorkerLoop',
    'UserError'
]


##################################################
## IMPORT
##################################################
def import_(*args):
    """Controller command that imports python modules on all workers.

    Each element of args should be a module name that is imported to
    all workers.

    Example:

    >>> pmi.import_('hello')
    >>> hw = pmi.create('hello.HelloWorld')
    """
    global inWorkerLoop
    if isController:
        if not args:
            raise UserError('pmi.import_ expects exactly 1 argument on controller!')

        # broadcast the statement
        _broadcast(_IMPORT, *args)
        # locally execute the statement
        return __workerImport_(*args)
    elif not inWorkerLoop:
        return receive(_IMPORT)


def __workerImport_(*modules):
    log.info("Importing modules: %s", modules)
    statement = 'import ' + ', '.join(modules)
    exec(statement, globals())


##################################################
## EXEC
##################################################
def exec_(*args):
    """Controller command that executes arbitrary python code on all workers.

    exec_() allows to execute arbitrary Python code on all workers.
    It can be used to define classes and functions on all workers.
    Modules should not be imported via exec_(), instead import_()
    should be used.

    Each element of args should be string that is executed on all
    workers.

    Example:

    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    """
    if __checkController(exec_):
        if not args:
            raise UserError('pmi.exec_ expects at least one argument(s) on controller!')

        # broadcast the statement
        _broadcast(_EXEC, *args)
        # locally execute the statement
        return __workerExec_(*args)
    else:
        return receive(_EXEC)


def __workerExec_(*statements):
    # executing the statement locally
    for statement in statements:
        log.info("Executing '%s'", statement)
        exec(statement, globals())


##################################################
## EXECFILE
##################################################
def execfile_(file):
    if __checkController(execfile_):
        _broadcast(_EXECFILE, file)
        return __workerExecfile_(file)
    else:
        return receive(_EXECFILE)


def __workerExecfile_(filename):
    log.info("Executing file '%s'", filename)
    exec(compile(open(filename, 'r').read(), filename, 'exec'), globals())


##################################################
## CREATE
##################################################
def create(cls=None, *args, **kwds):
    """Controller command that creates an object on all workers.

    cls describes the (new-style) class that should be instantiated.
    args are the arguments to the constructor of the class.  Only
    classes that are known to PMI can be used, that is, classes that
    have been imported to pmi via `exec_()` or `import_()`.

    Example:

    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    >>> print(hw)
    MPI process #0: Hello World!
    MPI process #1: Hello World!
    ...

    Alternative:
    Note that in this case the class has to be imported to the
    calling module *and* via PMI.

    >>> import hello
    >>> pmi.exec_('import hello')
    >>> hw = pmi.create(hello.HelloWorld)
    >>> print(hw)
    MPI process #0: Hello World!
    MPI process #1: Hello World!
    ...
    """
    if __checkController(create):
        if cls is None:
            raise UserError('pmi.create expects at least 1 argument on controller!')
        cls = _translateClass(cls)

        # generate a new oid
        oid = __OID()

        # translate the arguments
        cargs, ckwds, targs, tkwds = __translateArgs(args, kwds)
        # broadcast creation to the workers
        _broadcast(_CREATE, cls, oid, *targs, **tkwds)

        obj = __workerCreate(cls, oid, *cargs, **ckwds)

        # On the controller, store the oid in the instance
        obj.__pmioid = oid
        # Create the destroyer so that the instances on the workers
        # are destroyed
        obj.__pmidestroyer = __Destroyer(oid)

        return obj
    else:
        return receive(_CREATE)


def __workerCreate(cls, oid, *targs, **tkwds):
    # backtranslate the arguments
    args, kwds = __backtranslateOIDs(targs, tkwds)
    log.info('Creating: %s [%s]'
             % (__formatCall(cls.__name__, args, kwds), oid))
    # create the worker instance
    obj = cls(*args, **kwds)

    if isWorker:
        # store the new object
        if oid in OBJECT_CACHE:
            raise InternalError("Object [%s] is already in OBJECT_CACHE!" % oid)
        OBJECT_CACHE[oid] = obj
    return obj


##################################################
## CLONE
##################################################
# If the class is picklable, a living instance can be converted into a pmi object

##################################################
## CALL (INVOKE WITHOUT RESULT)
##################################################
def call(*args, **kwds):
    """Call a function on all workers, returning only the return value on the controller.

    function denotes the function that is to be called, args and kwds
    are the arguments to the function. If kwds contains keys that
    start with with the prefix '__pmictr_', they are stripped of the
    prefix and are passed only to the controller.
    If the function should return any results, it will be locally
    returned.
    Only functions that are known to PMI can be used, that is functions
    that have been imported to pmi via `exec_()` or `import_()`.

    Example:

    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    >>> pmi.call(hw.hello)
    >>> # equivalent:
    >>> pmi.call('hello.HelloWorld', hw)

    Note, that you can use only functions that are know to PMI when
    `call()` is called, i.e. functions in modules that have
    been imported via `exec_()`.
    """
    if __checkController(call):
        if len(args) == 0:
            raise UserError('pmi.call expects at least 1 argument on controller!')
        cfunction, tfunction, args = __translateFunctionArgs(*args)
        cargs, ckwds, targs, tkwds = __translateArgs(args, kwds)
        _broadcast(_CALL, tfunction, *targs, **tkwds)
        log.info("Calling: %s", __formatCall(cfunction, cargs, ckwds))
        return cfunction(*cargs, **ckwds)
    else:
        return receive(_CALL)


def __workerCall(function, *targs, **tkwds):
    function = __backtranslateFunctionArg(function)
    args, kwds = __backtranslateOIDs(targs, tkwds)
    log.info("Calling: %s", __formatCall(function, args, kwds))
    return function(*args, **kwds)


##################################################
## LOCAL CALL
##################################################
# provided for convenience
# you can just use: pmiobject.function(*args, **kwds)
def localcall(*args, **kwds):
    if __checkController(localcall):
        cfunction, tfunction, args = __translateFunctionArgs(*args)
        args, kwds = __translateProxies(args, kwds)
        log.info("Calling locally: %s", __formatCall(cfunction, args, kwds))
        return cfunction(*args, **kwds)
    else:
        raise UserError('Cannot call localcall on worker!')


##################################################
## INVOKE
##################################################
def invoke(*args, **kwds):
    """Invoke a function on all workers, gathering the return values into a list.

    function denotes the function that is to be called, args and
    kwds are the arguments to the function. If kwds contains keys
    that start with with the prefix '__pmictr_', they are stripped of
    the prefix and are passed only to the controller.

    On the controller, invoke() returns the results of the different
    workers as a list. On the workers, invoke returns None.
    Only functions that are known to PMI can be used, that is functions
    that have been imported to pmi via `exec_()` or `import_()`.

    Example:

    >>> pmi.exec_('import hello')
    >>> hw = pmi.create('hello.HelloWorld')
    >>> messages = pmi.invoke(hw.hello())
    >>> # alternative:
    >>> messages = pmi.invoke('hello.HelloWorld.hello', hw)
    """
    if __checkController(invoke):
        if len(args) == 0:
            raise UserError('pmi.invoke expects at least 1 argument on controller!')
        cfunction, tfunction, args = __translateFunctionArgs(*args)
        cargs, ckwds, targs, tkwds = __translateArgs(args, kwds)
        _broadcast(_INVOKE, tfunction, *targs, **tkwds)
        log.info("Invoking: %s", __formatCall(cfunction, cargs, ckwds))
        value = cfunction(*cargs, **ckwds)
        return _MPIGather(value)
    else:
        return receive(_INVOKE)


def __workerInvoke(function, *targs, **tkwds):
    function = __backtranslateFunctionArg(function)
    args, kwds = __backtranslateOIDs(targs, tkwds)
    log.info("Invoking: %s", __formatCall(function, args, kwds))
    value = function(*args, **kwds)
    return _MPIGather(value)


##################################################
## REDUCE (INVOKE WITH REDUCED RESULT)
##################################################
def reduce(*args, **kwds):
    """Invoke a function on all workers, reducing the return values to
    a single value.

    reduceOp is the (associative) operator that is used to process the
    return values, function denotes the function that is to be called,
    args and kwds are the arguments to the function. If kwds
    contains keys that start with with the prefix '__pmictr_', they
    are stripped of the prefix and are passed only to the controller.

    reduce() reduces the results of the different workers into a
    single value via the operation reduceOp. reduceOp is assumed to be
    associative.
    Both reduceOp and function have to be known to PMI, that is they
    must have been imported to pmi via `exec_()` or `import_()`.

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
    if __checkController(reduce):
        if len(args) <= 1:
            raise UserError('pmi.reduce expects at least 2 argument on controller!')
        # handle reduceOp argument
        creduceOp, treduceOp, args = __translateReduceOpArgs(*args)
        cfunction, tfunction, args = __translateFunctionArgs(*args)
        cargs, ckwds, targs, tkwds = __translateArgs(args, kwds)
        _broadcast(_REDUCE, treduceOp, tfunction, *targs, **tkwds)
        log.info("Reducing: %s", __formatCall(cfunction, cargs, ckwds))
        if _PMIComm and _PMIComm.isActive():
            if CONTROLLER in _PMIComm.getMPIcpugroup():
                value = cfunction(*args, **ckwds)
            else:
                value = None
        else:
            value = cfunction(*args, **ckwds)
        log.info("Reducing results via %s", creduceOp)
        return _MPIReduce(op=creduceOp, value=value)
    else:
        return receive(_REDUCE)


def __workerReduce(reduceOp, function, *targs, **tkwds):
    reduceOp = __backtranslateReduceOpArg(reduceOp)
    function = __backtranslateFunctionArg(function)
    args, kwds = __backtranslateOIDs(targs, tkwds)
    log.info("Reducing: %s", __formatCall(function, args, kwds))
    value = function(*args, **kwds)
    log.info("Reducing results via %s", reduceOp)
    return _MPIReduce(op=reduceOp, value=value)


##################################################
## SYNC
##################################################
def sync():
    """Controller command that deletes the PMI objects on the
    workers that have already been deleted on the controller.
    """
    if __checkController(sync):
        _broadcast(_SYNC)
    else:
        receive(_SYNC)


def __workerSync():
    """Worker sync is a nop, it only exists for the possible deletion
    of objects.
    """
    pass


##################################################
## DUMP
##################################################
def dump():
    """Controller function that dumps the object cache of PMI. For
    debugging purposes."""
    if __checkController(dump):
        _broadcast(_DUMP)
    else:
        receive(_DUMP)


def __workerDump():
    import pprint
    print(("OBJECT_CACHE=%s", pprint.pformat(OBJECT_CACHE)))


##################################################
## AUTOMATIC OBJECT DELETION
##################################################
def __delete():
    """Internal implementation of sync()."""
    global DELETED_OIDS
    if len(DELETED_OIDS) > 0:
        log.debug("Got %d objects in DELETED_OIDS.", len(DELETED_OIDS))
        __broadcastCmd(_DELETE, *DELETED_OIDS)
        DELETED_OIDS = []


def __workerDelete(*args):
    """Deletes the OBJECT_CACHE reference to a PMI object."""
    if len(args) > 0:
        log.info("Deleting oids: %s", args)
        for oid in args:
            obj = OBJECT_CACHE[oid]
            log.debug("  %s [%s]" % (obj, oid))
            # Delete the entry from the cache
            del OBJECT_CACHE[oid]


##################################################
## WORKER LOOP CONTROL
##################################################
def startWorkerLoop():
    """Worker command that starts the main worker loop.

    This function starts a loop that expects to receive PMI commands
    until `stopWorkerLoop()` or `finalizeWorkers()` is called on the
    controller.
    """
    global inWorkerLoop

    # On the controller, leave immediately
    if isController:
        log.info('Entering and leaving the worker loop')
        return None

    log.info('Entering the worker loop.')
    inWorkerLoop = True

    try:
        while True:
            receive()
    except StopIteration:
        inWorkerLoop = False


def finalizeWorkers():
    """Controller command that stops and exits all workers.
    """
    stopWorkerLoop(doExit=True)


def stopWorkerLoop(doExit=False):
    """Controller command that stops all workers.

    If doExit is set, the workers exit afterwards.
    """
    if __checkController(stopWorkerLoop):
        log.info('Calling all workers to stop.')
        _broadcast(_STOP, doExit)
    else:
        raise UserError('Cannot call stopWorkerLoop on worker!')


def __workerStop(doExit):
    if doExit:
        log.info('Stopping worker loop and exiting worker thread.')
        sys.exit()
    else:
        log.info('Stopping worker loop.')
        raise StopIteration()


def registerAtExit():
    """Controller command that registers the function
    `finalizeWorkers()` via atexit.
    """
    if __checkController(registerAtExit):
        import atexit
        atexit.register(finalizeWorkers)
    else:
        raise UserError('Cannot call registerAtExit on worker!')


##################################################
## PROXY METACLASS
##################################################

from types import BuiltinFunctionType, BuiltinMethodType, FunctionType, MethodType, LambdaType
from functools import partial


def is_function(obj):
    return isinstance(obj, (BuiltinFunctionType, BuiltinMethodType, FunctionType, MethodType, LambdaType, partial))


def _create_initializer(pmiobjectclassdef):
    def init(method_self, *args, **kwds):
        method_self.pmiobjectclassdef = pmiobjectclassdef
        pmiobjectclass = _translateClass(pmiobjectclassdef)
        method_self.pmiobject = create(pmiobjectclass, *args, **kwds)
        method_self.pmiobject._pmiproxy = method_self

    return init


def _create_caller(method_name):
    def caller(method_self, *args, **kwargs):
        method = getattr(method_self.pmiobject, method_name)
        return _backtranslateProxy(call(method, *args, **kwargs))

    return caller


def _create_invoker(method_name):
    def invoker(method_self, *args, **kwargs):
        method = getattr(method_self.pmiobject, method_name)
        return list(map(_backtranslateProxy, invoke(method, *args, **kwargs)))

    return invoker


def _create_local_caller(method_name):
    def local_caller(method_self, *args, **kwargs):
        method = getattr(method_self.pmiobject, method_name)
        return _backtranslateProxy(localcall(method, *args, **kwargs))

    return local_caller


def _create_property_local_getter(property_name):
    def property_getter(method_self):
        prop = getattr(method_self.pmiobject.__class__, property_name)
        fget = getattr(prop, 'fget')
        return _backtranslateProxy(fget(method_self.pmiobject))

    return property_getter


def _create_property_setter(property_name):
    def property_setter(method_self, value):
        setter_name = '.'.join((
            method_self.pmiobjectclassdef,
            property_name,
            'fset'
        ))
        return _backtranslateProxy(call(setter_name, method_self, value))

    return property_setter


class Proxy(type):
    """A metaclass to be used to create frontend serial objects."""

    def __init__(cls, name, bases, ns):
        if 'pmiproxydefs' in ns:
            defs = ns['pmiproxydefs']
            if 'cls' in defs:
                for base_class in bases:
                    if not hasattr(base_class, 'pmiproxydefs'):
                        continue
                    for k, v in base_class.pmiproxydefs.items():
                        if k == 'cls':
                            continue
                        if k in defs:
                            if isinstance(defs[k], Iterable):
                                defs[k] = list(defs[k]) + v
                        else:
                            defs[k] = v
                pmiobjectclassdef = defs['cls']
                proxy_init = _create_initializer(pmiobjectclassdef)
                setattr(cls, 'pmiinit', proxy_init)
                valid_init = isinstance(cls.__init__, (types.MethodType, types.FunctionType))
                if not valid_init:
                    cls.__init__ = proxy_init

            for method_name in defs.get('localcall', []):
                setattr(cls, method_name, _create_local_caller(method_name))

            for method_name in defs.get('pmicall', []):
                setattr(cls, method_name, _create_caller(method_name))

            for method_name in defs.get('pmiinvoke', []):
                setattr(cls, method_name, _create_invoker(method_name))

            for prop_name in defs.get('pmiproperty', []):
                newprop = property(
                    _create_property_local_getter(prop_name),
                    _create_property_setter(prop_name))
                setattr(cls, prop_name, newprop)


##################################################
## CONSTANTS AND EXCEPTIONS
##################################################

class InternalError(Exception):
    """Raised when PMI has encountered an internal error.

    Hopefully, this exceptions is never raised."""

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return workerStr + ': ' + self.msg

    def __repr__(self):
        return str(self)


class UserError(Exception):
    """Raised when PMI has encountered a user error.
    """

    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return workerStr + ': ' + self.msg

    def __repr__(self):
        return str(self)


##################################################
## BROADCAST AND RECEIVE
##################################################
def _broadcast(cmd, *args, **kwds):
    """Internal controller command that actually broadcasts a PMI command.

    The function first checks whether cmd is a valid PMI command, then
    it checks whether any objects have to be deleted before the
    command is broadcast, and finally it broadcasts the command
    itself.
    """
    __delete()
    log.debug("Broadcasting command: %s", _CMD[cmd][0])
    __broadcastCmd(cmd, *args, **kwds)


def __broadcastCmd(cmd, *args, **kwds):
    """This wraps a command with its argument into an internal __CMD
    object, so that it can be safely sent via MPI. __CMD is
    pciklable."""
    if not _checkCommand(cmd):
        raise InternalError('_broadcast needs a PMI command as first argument. Got %s instead!' % cmd)
    cmdobj = __CMD(cmd, args, kwds)
    _MPIBroadcast(cmdobj)


def receive(expected=None):
    """Worker command that receives and handles the next PMI command.

    This function waits to receive and handle a single PMI command. If
    expected is not None and the received command does not equal
    expected, raise a `UserError`.
    """
    __checkWorker(receive)
    if expected is None:
        log.debug('Waiting for next PMI command.')
    else:
        log.debug('Waiting for PMI command %s.', _CMD[expected][0])
    message = _MPIBroadcast()
    log.debug("Received message: %s", message)
    if type(message) is not __CMD:
        raise UserError("Received an MPI message that is not a PMI command: '%s'" % str(message))
    cmd = message.cmd
    args = message.args
    kwds = message.kwds
    if cmd == _DELETE:
        # if _DELETE is sent, delete the objects
        __workerDelete(*args)
        # recursively call receive once more
        return receive(expected)
    elif expected is not None and cmd != expected:
        # otherwise test whether the command is expected
        raise UserError("Received PMI command %s but expected %s" % (_CMD[cmd][0], _CMD[expected][0]))
    # determine which function to call
    cmd_func = _CMD[cmd][1]
    log.debug("Calling function %s", __formatCall(cmd_func.__name__, args, kwds))
    return cmd_func(*args, **kwds)


##################################################
## INTERNAL FUNTIONS
##################################################
class __OID:
    """Internal class that represents a PMI object id.

    An instance of this class can be pickled, so that it can be sent
    via MPI, and it is hashable, so that it can be used as a hash key
    (for OBJECT_CACHE).
    """

    def __init__(self):
        self.id = id(self)

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


class __Destroyer:
    def __init__(self, oid):
        self.oid = oid

    def __del__(self):
        log.info("Adding OID to DELETED_OIDS: [%s]", self.oid)
        DELETED_OIDS.append(self.oid)


class __CMD:
    """Internal, picklable class that represents a PMI
    command. Intended to be sent via MPI.
    """

    def __init__(self, cmd, args=None, kwds=None):
        if not _checkCommand(cmd):
            raise InternalError('Created __CMD object with invalid PMI command %s' % cmd)
        self.cmd = cmd
        self.args = args
        self.kwds = kwds

    def __str__(self):
        sargs = [_CMD[self.cmd][0]]
        if hasattr(self, 'args'):
            sargs.append(str(list(map(str, self.args))))
        if hasattr(self, 'kwds'):
            sargs.append(str(self.kwds))
        return 'PMICMD(%s)' % (', '.join(sargs))

    def __getstate__(self):
        state = (self.cmd, self.args, self.kwds)
        return state

    def __setstate__(self, state):
        self.cmd, self.args, self.kwds = state


def _isProxy(obj):
    return hasattr(obj, 'pmiobject')


def _checkCommand(cmd):
    return 0 <= cmd < _MAXCMD


def __checkController(func):
    """Checks whether we are on the controller, raises a UserError if
    we are on a worker and in the worker loop.

    Returns whether we are on the controller.
    """
    global inWorkerLoop
    if isController:
        return True
    else:
        if not inWorkerLoop:
            return False
        else:
            raise UserError("Cannot call %s on worker while in worker loop!" % func.__name__)


def __checkWorker(func):
    """Checks whether we are on a worker, raises a UserError if we are not.
    """
    if isController:
        raise UserError("Cannot call %s on the controller!" % func.__name__)


def _translateClass(cls):
    """Returns the class object of the class described by cls.
    """
    if cls is None:
        raise UserError("pmi.create expects at least 1 argument on controller")
    elif isinstance(cls, str):
        return eval(cls)
    elif isinstance(cls, type):
        return cls
    else:
        raise ValueError("__translateClass expects class as argument, but got %s" % cls)


def __mapArgs(func, args, kwds):
    """Internal function that maps a function to the list args and to
    the values of the dict kwds. Used by __translateArgs and
    __backtranslateOIDs.
    """
    targs = list(map(func, args))
    tkwds = {}
    for k, v in list(kwds.items()):
        tkwds[k] = func(v)
    return targs, tkwds


def _translateOID(obj):
    """Internal function that translates obj into an __OID
    object if it is a PMI object instance.

    If the object is not a PMI object, returns obj untouched.
    """
    if hasattr(obj, '__pmioid'):
        # if it is a pmi object, send the oid
        return obj.__pmioid
    else:
        return obj


def _backtranslateProxy(obj):
    # print('backTranslateProxy', type(obj))
    # Zobacz w starym espressopp jakie byly tu wartosci. Teraz ten pmi proxy zawiera pmi.Proxy typ
    # ktory przeciez jest metaklasa..
    if hasattr(obj, '_pmiproxy'):
        # print('pmiproxy type', type(obj._pmiproxy), dir(obj._pmiproxy))
        return obj._pmiproxy
    else:
        return obj


def _translateProxy(obj):
    if _isProxy(obj):
        return obj.pmiobject
    else:
        return obj


def __translateProxies(args, kwds):
    return __mapArgs(_translateProxy, args, kwds)


def __translateOIDs(args, kwds):
    """Internal function that translates all PMI object instances that
    occur in args or kwds into __OID objects that can be sent to the
    workers.
    """
    return __mapArgs(_translateOID, args, kwds)


def __translateArgs(args, kwds):
    args, kwds = __translateProxies(args, kwds)

    workerKwds = {}
    controllerKwds = {}
    for k in list(kwds.keys()):
        if k.startswith('__pmictr_'):
            knew = k[9:]
            controllerKwds[knew] = kwds[k]
        else:
            v = kwds[k]
            workerKwds[k] = v
            if k not in controllerKwds:
                controllerKwds[k] = v

    targs, tWorkerKwds = __translateOIDs(args, workerKwds)

    return args, controllerKwds, targs, tWorkerKwds


def _backtranslateOID(obj):
    """Internal worker function that backtranslates an __OID object
    into the corresponding PMI worker instance.

    If the object is not an __OID object, returns the object untouched.
    """
    if type(obj) is __OID:
        if obj in OBJECT_CACHE:
            return OBJECT_CACHE[obj]
        elif isController:
            # TODO: Not nice! Is this just so that broadcast can
            # return anything on the controller?
            return None
        else:
            raise InternalError("Object [%s] is not in OBJECT_CACHE" % obj)
    else:
        return obj


def __backtranslateOIDs(targs, tkwds):
    """Internal function that backtranslates all __OID object
    instances that occur in args or kwds into the cached PMI objects
    on the worker.
    """
    return __mapArgs(_backtranslateOID, targs, tkwds)


# Wrapper that allows to pickle a method
class __Method:
    def __init__(self, funcname, im_self, im_class=None):
        self.funcname = funcname
        self.im_self = _translateProxy(im_self)
        if im_class is None:
            self.im_class = self.im_self.__class__
        else:
            self.im_class = im_class
        self.__determineMethod()

    def __getstate__(self):
        return (self.funcname,
                _translateOID(self.im_self),
                self.im_class)

    def __setstate__(self, state):
        self.funcname, self.im_self, self.im_class = state
        self.im_self = _backtranslateOID(self.im_self)
        self.__determineMethod()

    def __determineMethod(self):
        for cls in self.im_class.mro():
            if hasattr(cls, self.funcname):
                function = getattr(cls, self.funcname)
                self.method = function.__get__(self.im_self, cls)
                break

    def __call__(self, *args, **kwds):
        return self.method(*args, **kwds)


# translate arguments to invoke
def __translateFunctionArgs(*args):
    """Internal controller function that normalizes the function
    argument to invoke(), call() or reduce().
    """
    if not args:
        raise TypeError("arguments missing")
    arg0 = args[0]
    if arg0 is None:
        raise TypeError("pmi expects function argument on controller")
    if isinstance(arg0, str):
        tfunction = arg0
        function = eval(arg0, globals())
        rargs = args[1:]
    elif isinstance(arg0, (types.FunctionType, types.BuiltinFunctionType)):
        if arg0.__name__ == '<lambda>':
            raise ValueError("pmi cannot handle lambda functions")
        tfunction = arg0
        function = tfunction
        rargs = args[1:]
    elif isinstance(arg0, (types.MethodType, types.BuiltinMethodType)):
        tfunction = __Method(arg0.__func__.__name__, arg0.__self__, arg0.__self__.__class__)
        function = tfunction
        rargs = args[1:]
    else:
        if len(args) <= 1:
            raise TypeError("got an object as first argument, but nothing as second")
        arg1 = args[1]
        if isinstance(arg1, str):
            tfunction = __Method(arg1, arg0)
            function = tfunction
            rargs = args[2:]
        else:
            raise ValueError("bad arguments")
    return function, tfunction, rargs


def __backtranslateFunctionArg(arg0):
    if isinstance(arg0, str):
        return eval(arg0, globals())
    else:
        return arg0


def __translateReduceOpArgs(*args):
    tfunction = _MPITranslateReduceOp(*args)
    if tfunction is not None:
        function = args[0]
        rargs = args[1:]
        return function, tfunction, rargs
    else:
        return __translateFunctionArgs(*args)


def __backtranslateReduceOpArg(arg0):
    function = _MPIBacktranslateReduceOp(arg0)
    if function is not None:
        return function
    else:
        return __backtranslateFunctionArg(arg0)


def __formatCall(function, args, kwds):
    def formatArgs(args, kwds):
        arglist = [repr(arg) for arg in args]
        for k, v in kwds.items():
            arglist.append('%s=%r' % (k, repr(v)))
        return ', '.join(arglist)

    return '%s(%s)' % (function, formatArgs(args, kwds))


# map of command names and associated worker functions
_CMD = [
    ('EXEC', __workerExec_),
    ('IMPORT', __workerImport_),
    ('EXECFILE', __workerExecfile_),
    ('CREATE', __workerCreate),
    ('INVOKE', __workerInvoke),
    ('CALL', __workerCall),
    ('REDUCE', __workerReduce),
    ('DELETE', __workerDelete),
    ('SYNC', __workerSync),
    ('STOP', __workerStop),
    ('DUMP', __workerDump)
]

_EXEC = 0
_IMPORT = 1
_EXECFILE = 2
_CREATE = 3
_INVOKE = 4
_CALL = 5
_REDUCE = 6
_DELETE = 7
_SYNC = 8
_STOP = 9
_DUMP = 10

_MAXCMD = len(_CMD)

# set that stores which oids have been deleted
DELETED_OIDS = []
# dict that stores the objects corresponding to an oid
OBJECT_CACHE = {}

inWorkerLoop = False

##################################################
## MPI SETUP
##################################################
from mpi4py import MPI
from mpi4py.MPI import OP_NULL, MAX, MIN, SUM, PROD, LAND, BAND, LOR, BOR, LXOR, BXOR, MAXLOC, MINLOC, REPLACE

log = logging.getLogger('%s.controller' % __name__)


# log.setLevel(logging.DEBUG)

def _MPIInit(comm=MPI.COMM_WORLD):
    # The communicator used by PMI
    global _MPIcomm, _PMIComm, CONTROLLER, rank, size, \
        isController, isWorker, workerStr, log

    _MPIcomm = comm
    _PMIComm = None

    CONTROLLER = 0
    rank = _MPIcomm.rank
    size = _MPIcomm.size

    # whether this is a worker or a controller
    isController = rank == CONTROLLER
    isWorker = not isController

    if isController:
        workerStr = 'Controller'
        log = logging.getLogger('%s.controller' % __name__)
    else:
        workerStr = 'Worker %d' % rank
        log = logging.getLogger('%s.worker%d' % (__name__, rank))


def _MPIGather(value):
    global CONTROLLER, _MPIcomm, _PMIComm
    if _PMIComm and _PMIComm.isActive():
        return _PMIComm.getMPIsubcommWithController().gather(value, root=CONTROLLER)
    else:
        return _MPIcomm.gather(value, root=CONTROLLER)


def _MPIBroadcast(value=None):
    global CONTROLLER, _MPIcomm, _PMIComm
    if _PMIComm and _PMIComm.isActive():
        return _PMIComm.getMPIsubcommWithController().bcast(value, root=CONTROLLER)
    else:
        return _MPIcomm.bcast(value, root=CONTROLLER)


def _MPIReduce(op, value):
    global CONTROLLER, _MPIcomm, _PMIComm
    if _PMIComm and _PMIComm.isActive():
        return _PMIComm.getMPIsubcommWithController().reduce(value, root=CONTROLLER, op=op)
    else:
        return _MPIcomm.reduce(value, root=CONTROLLER, op=op)


def _MPISpawnAndMerge(ntasks, command):
    if _PMIComm and _PMIComm.isActive():
        raise UserError('Requested to SpawnAndMerge, but worker subcommunicator group is active.')

    cmd = command[0]
    if len(command) > 1:
        args = command[1:]
    else:
        args = ()
    intercomm = \
        _MPIcomm.Spawn(cmd, args,
                       maxprocs=ntasks, root=CONTROLLER)
    newcomm = intercomm.Merge(False)
    _MPIInit(newcomm)


def _MPIMergeWithParent():
    if _PMIComm and _PMIComm.isActive():
        raise UserError('Requested to MergeWithParent, but worker subcommunicator group is active.')

    intercomm = _MPIComm.Get_parent()
    if intercomm == MPI.COMM_NULL: return
    newcomm = intercomm.Merge(True)
    _MPIInit(newcomm)


# map of command names and associated worker functions

_REDUCEOP = [OP_NULL, MAX, MIN, SUM, PROD, LAND, BAND, LOR, BOR,
             LXOR, BXOR, MAXLOC, MINLOC, REPLACE]


class _ReduceOp(object):
    def __init__(self, op):
        self.op = op

    def __getstate__(self):
        i = 0
        for op in _REDUCEOP:
            if self.op is op:
                return i
            i += 1

    def __setstate__(self, state):
        self.op = _REDUCEOP[state]

    def getOp(self):
        return self.op


def _MPITranslateReduceOp(*args):
    arg0 = args[0]
    if isinstance(arg0, MPI.Op):
        return _ReduceOp(arg0)
    else:
        return None


def _MPIBacktranslateReduceOp(arg0):
    if isinstance(arg0, _ReduceOp):
        return arg0.getOp()
    else:
        return None


########################################
# Communicator Class
########################################

class CommunicatorLocal(object):
    'PMI CommunicatorLocal class'

    def __init__(self, cpugroup=None):
        if not cpugroup:
            cpugroup = range(0, _MPIcomm.size)
        self._cpugroup = cpugroup
        self._MPIsubcomm = None
        self._MPIsubcommWithController = None
        self._isActive = False
        if (max(self._cpugroup) < _MPIcomm.size) and (min(self._cpugroup) >= 0):
            commgroup = _MPIcomm.Get_group()
            subcommgroup = commgroup.Incl(self._cpugroup)
            self._MPIsubcomm = _MPIcomm.Create(subcommgroup)
            if CONTROLLER in self._cpugroup:
                self._MPIsubcommWithController = self._MPIsubcomm
            else:
                cpugroupWithController = self._cpugroup[:]
                cpugroupWithController.insert(0, CONTROLLER)
                commgroupWithController = _MPIcomm.Get_group()
                subcommgroupWithController = commgroupWithController.Incl(cpugroupWithController)
                self._MPIsubcommWithController = _MPIcomm.Create(subcommgroupWithController)

    def getMPIsubcomm(self):
        'getter for MPIsubcomm'
        return self._MPIsubcomm

    def getMPIsubcommWithController(self):
        'getter for MPIsubcomm'
        return self._MPIsubcommWithController

    def getMPIcpugroup(self):
        'getter for MPIsubgroup'
        return self._cpugroup

    def isActive(self):
        return self._isActive

    def activate(self):
        self._isActive = True

    def deactivate(self):
        self._isActive = False


class Communicator(object):
    """PMI Communicator class"""

    def __init__(self, *args, **kwds):
        if isController:
            self.localcomm = create(_translateClass(CommunicatorLocal), *args, **kwds)

    def getMPIsubcomm(self):
        return self.localcomm._MPIsubcomm

    def getMPIsubcommWithController(self):
        return self.localcomm._MPIsubcommWithController

    def getMPIcpugroup(self):
        return self.localcomm._cpugroup

    def isActive(self):
        return self.localcomm._isActive

    def activate(self):
        self.localcomm._isActive = True

    def deactivate(self):
        self.localcomm._isActive = False


##################################################
## check if worker is active
##################################################
def workerIsActive():
    if _PMIComm and _PMIComm.isActive():
        if _MPIcomm.rank in _PMIComm.getMPIcpugroup():
            return True
        else:
            return False
    else:
        return True


##################################################
## SETUP
##################################################
def setup(ntasks=None,
          taskcmd=(sys.argv[0], os.path.abspath('pmi.py'))):
    if isController:
        registerAtExit()
        if ntasks is not None:
            if size > 1:
                # parallel tasks already running
                if size != ntasks:
                    raise UserError(
                        'setup() requested to start %d tasks, but %d tasks are already running.' % (ntasks, size))
            else:
                # start ntasks tasks
                _MPISpawnAndMerge(ntasks - 1, taskcmd)

    else:
        startWorkerLoop()


##################################################
## MODULE BODY
##################################################
# set up MPI
_MPIInit()

# if the module is executed as script, it was probably spawned
if __name__ == 'main':
    _MPIMergeWithParent()
    startWorkerLoop()
