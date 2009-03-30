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
    def __init__(self, arg=None) :
        object.__init__(self)
        self.arg = arg

    def f(self) :
        self.f_arg = 42
        return 42

    def g(self, a) :
        self.g_arg = a
        return a
