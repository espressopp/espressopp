import espresso
import sandbox

class MyA(sandbox.A):
    def f(self):
        print("I'm PyA.f()")

class MyBsimple(sandbox.Bsimple):
    def f(self):
        print("I'm PyBsimple.f()")

class MyB(sandbox.B):
    def f(self):
        print("I'm PyB.f()")

a = MyA()
a.g()

bsimple=MyBsimple()
bsimple.g()

b=MyB()
b.g()

sandbox.caller(b)
