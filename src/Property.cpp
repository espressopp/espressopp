#include "python.hpp"
#include "Property.hpp"
#include "Particle.hpp"

using namespace espresso;
using namespace espresso::particles;

PropertyBase::
PropertyBase(const Storage::SelfPtr _storage)
  : storage(_storage) {}

Storage::SelfPtr
const PropertyBase::getStorage() const { return storage; }

void espresso::registerPythonProperty()
{
  using namespace espresso::python;
  class_< PropertyBase, boost::noncopyable >
    ("PropertyBase", no_init)
    ;

  class_< Property< Real3D >, bases< PropertyBase > >
    ("Real3DProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< Real3D >::getItem)
    .def("__setitem__", &Property< Real3D >::setItem);
  
  class_< Property< int >, bases< PropertyBase > >
    ("IntegerProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< int >::getItem)
    .def("__setitem__", &Property< int >::setItem);
  
  class_< Property< real >, bases< PropertyBase > >
    ("RealProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< real >::getItem)
    .def("__setitem__", &Property< real >::setItem);
  
  class_< ArrayProperty< int >, bases< PropertyBase > >
    ("IntegerArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &ArrayProperty< int >::getItem)
    .def("__setitem__", &ArrayProperty< int >::setItem);
  
  class_< ArrayProperty< real >, bases< PropertyBase > >
    ("RealArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &ArrayProperty< real >::getItem)
    .def("__setitem__", &ArrayProperty< real >::setItem);
}

