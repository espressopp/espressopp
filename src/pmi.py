import types
import logging
from espresso import esmpi as mpi

log = logging.getLogger('espresso.pmi')

def exec_(statement) :
    """This allows to import modules and to define classes and
    functions on all workers.
    """
    log.info("executing '%s'", statement)
    exec statement in globals()

# import_ is an alias for exec_
import_=exec_

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
        # the class is passed as name
        pass
    elif type(theClass) == types.TypeType :
        # a class object is passed
        theClass = theClass.__module__ + '.' + theClass.__name__
    elif type(theClass) == types.ClassType :
        raise TypeError("pmi can't handle old-style classes. \
Please create old style classes via their names.")
    else :
        raise TypeError("expects class as first argument")

    # create the class creation string
    s = "%s%s" % (theClass, args)
    log.info("creating class '%s'", s)
    # broadcast it to all workers
    return eval(s)

def invoke(function, *args) :
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
    log.info("invoking %s", function.__name__)
    return function(*args)

def delete(obj) :
    del(obj)

def finish() :
    'Ends PMI, i.e. stops all workers'
    pass

def workerLoop() :
    'Starts the main worker loop. Has to be called on the workers.'
    pass

