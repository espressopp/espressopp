# Fake PMI
# This just shows what interface I expect that can be used

def importModule(name) :
    return __import__(name)

def create(theClass, *args) :
    """Create an object on all workers.

    theClass is the class that should be used, *args are the arguments
    to the constructor of the class.

    Note, that you can use only those classes that are know to PMI
    when this function is called, i.e. classes in modules that have
    been imported via pmi.importModule().
    """
    return theClass(*args)

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
    return function(*args)

def eval(code) :
    """Even this might be possible! This should allow to define Python
    classes and functions on all workers.
    """
    pass

def delete(obj) :
    del(obj)

def finish() :
    'Ends PMI, i.e. stops all workers'
    pass

def _workerLoop() :
    'Starts the main worker loop. Has to be called on the workers.'
    pass

# Call the workerLoop
_workerLoop()
