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

from _espresso import esutil_Real3D 
class Real3D (object) :
    # TODO: make flexible: allow to be initialized from any sequence
    def __init__(self, x=0.0, y=0.0, z=0.0) :
        object.__init__(self)
        self.__cc = esutil_Real3D(x, y, z)

    # wrap setitem and getitem
    def __getitem__(self, i) :
        return self.__cc[i]
    def __setitem__(self, i, v) :
        self.__cc[i] = v

    # create setters and getters
    @propget
    def x(self) :
        return self.__cc[0]

    @propset
    def x(self, v) :
        self.__cc[0] = v

    @propget
    def y(self) :
        return self.__cc[1]

    @propset
    def y(self, v) :
        self.__cc[1] = v

    @propget
    def z(self) :
        return self.__cc[2]

    @propset
    def z(self, v) :
        self.__cc[2] = v

    # string conversion
    def __str__(self) :
        return str(tuple(self))

    def __add__(self, o) :
        return self.__cc + o.__cc

    def __mul__(self, o) :
        return self.__cc * o.__cc
