#include "Real3D.hpp"
#include <boost/python.hpp>

void 
espresso::esutil::registerPythonReal3D() {
  using namespace espresso::esutil;
  using namespace boost::python;

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
