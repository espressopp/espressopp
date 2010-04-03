import poc
from _call_policies import *

# When creating an A, it is kept in a shared_ptr
print('* A as shared_ptr *')
a = A()
print(a.s)
b = a
a.s = "Hello, Olaf!"
print(b.s)
del(a)
print(b.s)

# Putting a into container
print('* Put A into container *')
a = A()
c = Container()
c.setA(a);
a2 = c.getA();
print(a2.s)
c.setA(a2)
del(c)

# When returning an A from a Container, 
# I can get it as internal reference
print('* A as internal reference *')
c = Container()
a = c.getA()
print(a.s)
del(c)
print(a.s)
del(a)

# When getting the A as an existing object,
# it will be destroyed when the Container is deleted
print('* A as reference to existing object *')
c = Container()
a_danger = c.getADanger()
print(a_danger.s)
del(c)
print(a_danger.s) # Segfault!

