import poc
from _diamond import A as _A, B as _B

class A(_A):
    # Call the C++ inherited constructor and set cxxinit.
    def __init__(self, _i):
        print('Python A constructor start')
        if not hasattr(self, 'cxxinit'):
            _A.__init__(self)
            self.cxxinit = True
        print('Python A constructor end')

    # Only the Python class where a certain C++ method is provided for
    # the first time calls the C++ method.
    def doIt(self, i):
        print('Python A doIt start')
        _A.doIt(self)
        print('Python A doIt end')

class B(A, _B):
    # First call the C++ inherited constructor and set cxxinit, only
    # then call the Python inherited constructor.
    def __init__(self):
        print('Python B constructor start')
        if not hasattr(self, 'cxxinit'):
            _B.__init__(self)
            self.cxxinit = True
        A.__init__(self, 52)
        print('Python B constructor end')

    # Do NOT call the C++ method, but the Python inherited method.
    def doIt(self):
        print('Python B doIt start')
        A.doIt(self, 42)
        print('Python B doIt end')

b = B()
b.doIt()
b.doIt2()
