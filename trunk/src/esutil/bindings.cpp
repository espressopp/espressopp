#include "bindings.hpp"
#include "Real3D.hpp"
#include <boost/python.hpp>

namespace espresso {
  namespace esutil {
    
    using namespace boost::python;
    /** Register all C++ classes of namespace esutil 
	
     * Real3D
     
     */
    void registerPython() {
      // Export the C++ class Real3D to Python 
      class_<Real3D>("esutil_Real3D", init<real, real, real>())
	.def("__getitem__", &Real3D::getItem)
	.def("__setitem__", &Real3D::setItem)
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
