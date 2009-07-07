#include "Real3D.hpp"
#include "python.hpp"

namespace espresso { 
  struct real3D_pickle_suite 
    : boost::python::pickle_suite
  {
    static
    espresso::python::tuple
    getinitargs(const Real3D& v)
    { return espresso::python::make_tuple(v[0], v[1], v[2]); }
  };
    
  void 
  registerPythonReal3D() {
    // here we do NOT use espresso::python
    using namespace boost::python;
      
    // Export the C++ class Real3D to Python 
    class_<Real3D>("Real3D", init<real, real, real>())
      .def("__getitem__", &Real3D::getItem)
      .def("__setitem__", &Real3D::setItem)
      .def(self * real())
      .def(real() * self)
      .def(self * self)
      .def(self + self)
      .def(self - self)
      .def(self += self)
      .def(self -= self)
      .def(self == self)
      .def(self != self)
      .def("cross", &Real3D::cross)
      .def("sqr", &Real3D::sqr)
      .def("__abs__", &Real3D::abs)
      .def_pickle(real3D_pickle_suite())
      ;
  }
}
