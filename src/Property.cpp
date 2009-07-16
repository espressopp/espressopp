#include "python.hpp"
#include "Property.hpp"
#include "Particle.hpp"

void espresso::registerPythonProperties()
{
  using namespace espresso::particles;
  using namespace espresso::python;

  class_< Property< Real3D > >
    ("Real3DProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< Real3D >::getItem)
    .def("__setitem__", &Property< Real3D >::setItem);
  
  class_< Property< int > >
    ("IntegerProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< int >::getItem)
    .def("__setitem__", &Property< int >::setItem);
  
  class_< Property< real > >
    ("RealProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< real >::getItem)
    .def("__setitem__", &Property< real >::setItem);
  
  class_< ArrayProperty< int > >
    ("IntegerArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &ArrayProperty< int >::getItem)
    .def("__setitem__", &ArrayProperty< int >::setItem);
  
  class_< ArrayProperty< real > >
    ("RealArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &ArrayProperty< real >::getItem)
    .def("__setitem__", &ArrayProperty< real >::setItem);
}

