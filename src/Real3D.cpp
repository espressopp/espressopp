#include "Real3D.hpp"
#include "python.hpp"

using namespace espresso::python;

namespace espresso { 
  struct real3D_pickle_suite : pickle_suite
  {
    static
    espresso::python::tuple
    getinitargs(const Real3D v)
    { return make_tuple(v[0], v[1], v[2]); }
  };
    
  void 
  registerPythonReal3D() {
    // here we do NOT use espresso::python::class_
    // Export the C++ class Real3D to Python 
    boost::python::class_< Real3D >
      ("Real3D", init<real, real, real>())
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
