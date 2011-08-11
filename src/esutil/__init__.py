"""This module defines helper functions for python.
Modify with care, it contains some **black magic**.

ExtendBaseClass
---------------

Use this class as a meta class for an object, and the base class of
the object will be extended by the functions of this class.
    
Example:

>>> class Test :
>>>   def test(): return 'Test'
>>>
>>> class TestExtender(Test) :
>>>   __metaclass__ = ExtendBaseClass
>>>   def test_extend(): return 'Test extension'
>>>
>>> t = Test()
>>> t.test()
>>> t.test_extend()
'Test'
'Test extension'

Stolen and modified from
http://code.activestate.com/recipes/412717/ and
http://www.boost.org/doc/libs/1_35_0/libs/python/doc/tutorial/doc/html/python/techniques.html#python.extending_wrapped_objects_in_python
"""
from espresso import pmi
if pmi.isController :
    def pmiimport(module):
        pmi.exec_('import ' + module)
else:
    def pmiimport(module):
        pass
        
pmiimport('espresso.esutil')

from espresso.esutil.RNG import *
from espresso.esutil.UniformOnSphere import *
from espresso.esutil.NormalVariate import *

class ExtendBaseClass (type) :
    def __new__(self, name, bases, dict):
        del dict['__metaclass__']
        del dict['__module__']

        theClass = bases[0]
        # loop over all items in the class and replace it
        for k,v in dict.iteritems():
            setattr(theClass, k, v)
        return theClass

def choose(val, altval) :
    if (val is None) :
        return altval
    else :
        return val

def cxxinit(obj, cls, *args, **kwds):
#     # check whether the class is a boost.python class
#     if not issubclass(cls, type):
#         raise TypeError('cxxinit requires a class input argument.')
    if not hasattr(obj, 'cxxclass'):
        obj.cxxclass = cls
        cls.__init__(obj, *args, **kwds)

# def pmiinit(obj, cls, *args, **kwds):
#     if not hasattr(obj, 'pmiobject'):
#         obj.pmiobject = pmi.create(cls, *args, **kwds)

