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

# Make the property setter decorator syntax of python 2.6+ available
# to earlier versions
try :
    __setter = property.setter
except AttributeError :
    import __builtin__, sys
    # save the property builtin
    _property = __builtin__.property
    # now define our property
    # stolen from http://bruynooghe.blogspot.com/2008/04/xsetter-syntax-in-python-25.html 
    class property(_property):
        def __init__(self, fget, *args, **kwargs):
            self.__doc__ = fget.__doc__
            super(property, self).__init__(fget, *args, **kwargs)

        def setter(self, fset):
            cls_ns = sys._getframe(1).f_locals
            for k, v in cls_ns.iteritems():
                if v == self:
                    propname = k
                    break
            cls_ns[propname] = _property(self.fget, fset,
                                        self.fdel, self.__doc__)
            return cls_ns[propname]

        def deleter(self, fdel):
            cls_ns = sys._getframe(1).f_locals
            for k, v in cls_ns.iteritems():
                if v == self:
                    propname = k
                    break
            cls_ns[propname] = _property(self.fget, self.fset,
                                        fdel, self.__doc__)
            return cls_ns[propname]

    # Now override the property builtin
    __builtin__.property = property
