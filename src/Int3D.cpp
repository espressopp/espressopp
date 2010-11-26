#include "Int3D.hpp"
#include <python.hpp>
#include <boost/python/implicit.hpp>

namespace espresso {
  struct int3D_pickle_suite : boost::python::pickle_suite {
    static
    python::tuple
    getinitargs(const Int3D v)
    { return python::make_tuple(v[0], v[1], v[2]); }
  };
  
  void
  Int3D::
  registerPython() {
    using namespace python;
    using namespace boost::python;
    class_< Int3D >("Int3D", init<>())
      .def(init< double, double, double >())
      .def("__getitem__", &Int3D::getItem)
      .def("__setitem__", &Int3D::setItem)
      .def(self += self)
      .def(self -= self)
      .def(self == self)
      .def(self != self)
      .def(self + self)
      .def(self - self)
      .def_pickle(int3D_pickle_suite())
      ;

  }
}
