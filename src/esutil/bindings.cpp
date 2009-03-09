#include "bindings.hpp"
#include <boost/python.hpp>
#include "Real3D.hpp"

namespace espresso {
  namespace esutil {
    
  using namespace boost::python;
    
    /** Help routine to translate a Real3D to a Python tuple with three doubles.
	
	\param x is a 3D real vector
	\returns a Python tuple with the three vector elements
    */
    
    object makeTupleReal3D(Real3D x) {
      return make_tuple(x[0], x[1], x[2]);
    }
    
    /** Register all C++ classes of namespace esutil 
	
     * Real3D
     
     */
    
    void registerPython() {
      
      // Export the C++ class Real3D to Python via Boost-Python
      
      class_<Real3D>("esutil_Real3D", init<real, real, real>())
	.def("__getitem__", &Real3D::getItem)
	.def("__setitem__", &Real3D::setItem)
	.def("tuple", &makeTupleReal3D)
	.def(self * real())
	.def(real() * self)
	.def(self * self)
	.def(self + self)
	.def(self - self)
	.def(self += self)
	.def(self -= self)
	.def("sqr", &Real3D::sqr)
	;
      
      // ToDo: __eq__, __ne__ ?
      
    }
  }
}
