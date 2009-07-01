import poc
from _diamond import A as _A, B as _B

class A(_A):
    def __init__(self, _i):
        print('Python A constructor start')
        _A.__init__(self)
        print('Python A constructor end')

    def doIt(self):
        super(A, self).doIt()
        print('Python A')

class B(A, _B):
    def __init__(self):
        print('Python B constructor start')
#        A.__init__(self, 1)
        _B.__init__(self)
        print('Python B constructor end')

    def doIt(self):
        super(B, self).doIt()
        print('Python B')

b = B()
b.doIt()
b.doIt2()
