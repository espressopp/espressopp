import poc
import _inheritance

class MyA(_inheritance.A):
    def f(self):
        print("I'm PyA.f()")

class MyBsimple(_inheritance.Bsimple):
    def f(self):
        print("I'm PyBsimple.f()")

class MyB(_inheritance.B):
    def f(self):
        print("I'm PyB.f()")

a = MyA()
a.g()

bsimple=MyBsimple()
bsimple.g()

b=MyB()
b.g()

_inheritance.caller(b)
