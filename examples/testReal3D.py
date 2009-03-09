#########################################
#
# ToDo: JH will make y Python UnitTest of this
#
#  + implement __eq__  - operator ==
#  + implement __ne__  - operator !=
#  + export -= to Python
#
#  + allow to pickle Real3D values to transmit with boost.mpi

print 'This is a test of Real3D'

from _espresso import Real3D

x = Real3D(1.2, 1.3, 1.4)

print 'x = ', x
print 'x = ', x.tuple()

print 'x = ', x[0], x[1], x[2]

try:
   print 'x[-1] = ', x[-1]
except IndexError:
   print 'IndexError is fine'

try:
   x[3] = 3.3
except IndexError:
   print 'IndexError is fine'

x[1] = 0.0
print 'x = ', x[0], x[1], x[2]
x = x * 3
print 'x = ', x[0], x[1], x[2]
x = 0.3 * x
print 'x = ', x[0], x[1], x[2]
dot = x * x
print 'dot x = ', dot

y = x  +  x
print 'x + x = ', y[0], y[1], y[2]
y = y  -  x
print 'y - x = ', y[0], y[1], y[2]

print 'x = ', x[0], x[1], x[2]
x += x
print 'x += x', x[0], x[1], x[2]

print 'sqr(x) = ', x.sqr()

