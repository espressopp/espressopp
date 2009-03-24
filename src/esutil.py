"""
This module defines helper functions for python.

It contains some **black magic**. Use this class as a meta class for an object,
and the base class of the object will be extended by the functions
of this class.

Example:

>>>class Test :
>>>  def test(): return 'Test'

>>>class TestExtender(Test) :
>>>  __metaclass__ = ExtendBaseClass
>>>  def test_extend(): return 'Test extension'

>>>t = Test()
>>>assert t.test() == 'Test'
>>>assert t.test_extend() == 'Test extension'

Stolen and modified from
http://code.activestate.com/recipes/412717/ and
http://www.boost.org/doc/libs/1_35_0/libs/python/doc/tutorial/doc/html/python/techniques.html#python.extending_wrapped_objects_in_python
"""

import sys

def propget(func):
    'Decorator to easily define property getters.'
    locals = sys._getframe(1).f_locals
    name = func.__name__
    prop = locals.get(name)
    if not isinstance(prop, property):
        prop = property(func, doc=func.__doc__)
    else:
        doc = prop.__doc__ or func.__doc__
        prop = property(func, prop.fset, prop.fdel, doc)
    return prop

def propset(func):
    'Decorator to easily define property setters.'
    locals = sys._getframe(1).f_locals
    name = func.__name__
    prop = locals.get(name)
    if not isinstance(prop, property):
        prop = property(None, func, doc=func.__doc__)
    else:
        doc = prop.__doc__ or func.__doc__
        prop = property(prop.fget, func, prop.fdel, doc)
    return prop

def propdel(func):
    'Decorator to easily define property deletion.'
    locals = sys._getframe(1).f_locals
    name = func.__name__
    prop = locals.get(name)
    if not isinstance(prop, property):
        prop = property(None, None, func, doc=func.__doc__)
    else:
        prop = property(prop.fget, prop.fset, func, prop.__doc__)
    return prop

class ExtendBaseClass(type):
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
