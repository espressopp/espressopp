#include "RealTensor.hpp"
#include "python.hpp"

using namespace espresso::python;

namespace espresso { 
  struct realTensor_pickle_suite : pickle_suite
  {
    static
    espresso::python::tuple
    getinitargs(const RealTensor v)
    { return make_tuple(v[0], v[1], v[2], v[3], v[4], v[5]); }
  };
    
  void 
  registerPythonRealTensor() {
    // here we do NOT use espresso::python::class_
    // Export the C++ class RealTensor to Python 
    boost::python::class_< RealTensor >
      ("RealTensor", init<real, real, real, real, real, real>())
      .def("__getitem__", &RealTensor::getItem)
      .def("__setitem__", &RealTensor::setItem)
      .def(self * real())
      .def(real() * self)
      .def(self * self)
      .def(self + self)
      .def(self - self)
      .def(self += self)
      .def(self -= self)
      .def(self == self)
      .def(self != self)
      .def("__abs__", &RealTensor::abs)
      .def_pickle(realTensor_pickle_suite())
      ;
  }
}
