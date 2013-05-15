#include <python.hpp>
#include <boost/python/implicit.hpp>
#include "Tensor.hpp"

namespace espresso {
  struct tensor_pickle_suite : boost::python::pickle_suite {
    static
    python::tuple
    getinitargs(const Tensor v)
    { return python::make_tuple(v[0], v[1], v[2], v[3], v[4], v[5]); }
  };
  
  void
  Tensor::
  registerPython() {
    using namespace python;
    using namespace boost::python;
    class_< Tensor >("Tensor", init<>())
      .def(init< double, double, double, double, double, double >())
      .def("__getitem__", &Tensor::getItem)
      .def("__setitem__", &Tensor::setItem)
      .def(self += self)
      .def(self -= self)
      .def(self *= real())
      .def(self /= real())
      .def(self == self)
      .def(self != self)
      .def(self + self)
      .def(self - self)
      .def(self * real())
      .def(self / real())
      .def(real() * self)
      .def_pickle(tensor_pickle_suite())
      ;

  }
}
