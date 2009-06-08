# Test module to test PMI imports

class AOS :
    pass

def f() :
    global f_arg
    f_arg = 42
    return 42

def g(a) :
    global g_arg
    g_arg = a
    return a

def add(a, b) :
    return a + b

class A(object) :
    created = False
    initarg = None
    initkwds = None
    f_called = False
    g_arg = None
    g_kwds = None
    
    def __init__(self, arg=None, **kwds) :
        A.created = True
        A.initarg = arg
        A.initkwds = kwds
        A.f_called = False
        A.g_arg = None
        A.g_kwds = None

    def f(self) :
        A.f_called = True
        return 42

    def g(self, arg=None, **kwds) :
        A.g_arg = arg
        A.g_kwds = kwds
        return arg

    def __del__(self):
        A.created = False
        A.initarg = None
        A.initkwds = None
        A.f_called = False
        A.g_arg = None
        A.g_kwds = None


class Local(object):
    called = None
    def __init__(self):
        Local.called = 'init'
    def f(self):
        Local.called = 'f'
        return 'f'
    def g(self):
        Local.called = 'g'
        return 'g'
    def h(self):
        Local.called = 'h'
        return 'h'
    @property
    def x(self):
        Local.called = 'x.get'
        return self._x
    @x.setter
    def x(self, val):
        Local.called = 'x.set'
        self._x = val
    def __del__(self):
        Local.called = None
