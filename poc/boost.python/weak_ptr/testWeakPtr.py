import poc
import _weak_ptr

a = _weak_ptr.A(3)
_weak_ptr.setWeakPtr(a)
_weak_ptr.outputPtr()
