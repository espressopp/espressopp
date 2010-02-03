#include "Int3D.hpp"
#include "Int3DRef.hpp"
#include <python.hpp>
#include <boost/python/implicit.hpp>

namespace espresso {
  Int3D::Int3D(const Int3DRef &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }

  Int3D::Int3D(const ConstInt3DRef &v) {
    for (int i = 0; i < 3; i++)
      data[i] = v[i];
  }

  Int3D::operator Int3DRef() {
    return Int3DRef(data);
  }
  
  Int3D::operator ConstInt3DRef() const {
    return ConstInt3DRef(data);
  }
  
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

    boost::python::implicitly_convertible<Int3D, Int3DRef>();
    boost::python::implicitly_convertible<Int3DRef, Int3D>();
    boost::python::implicitly_convertible<const Int3D, ConstInt3DRef>();
    boost::python::implicitly_convertible<ConstInt3DRef, Int3D>();

  }
}
