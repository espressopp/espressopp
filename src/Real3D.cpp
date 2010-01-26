#include "Real3D.hpp"
#include "types.hpp"
#include <python.hpp>
#include <boost/python/implicit.hpp>

namespace espresso {
  Real3D::Real3D(const Real3DRef &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }

  struct real3D_pickle_suite : boost::python::pickle_suite {
    static
    python::tuple
    getinitargs(const Real3D v)
    { return python::make_tuple(v[0], v[1], v[2]); }
  };
  
  void
  Real3D::
  registerPython() {
    using namespace python;
    using namespace boost::python;
    class_< Real3D >("Real3D", init<>())
      .def(init< double, double, double >())
      .def("__getitem__", &Real3D::getItem)
      .def("__setitem__", &Real3D::setItem)
      .def("sqr", &Real3D::sqr)
      .def("__abs__", &Real3D::abs)
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
      .def(self * self)
      .def("cross", &Real3D::cross)
      .def_pickle(real3D_pickle_suite())
      ;

    boost::python::implicitly_convertible<Real3D, Real3DRef>();
    boost::python::implicitly_convertible<Real3DRef, Real3D>();

  }
}
