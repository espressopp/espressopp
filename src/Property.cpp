#include "python.hpp"
#include "Property.hpp"
#include "Particle.hpp"

using namespace espresso;
using namespace espresso::particles;
using namespace espresso::storage;

void espresso::registerPythonProperty()
{
  using namespace espresso::python;
  class_< PropertyBase, boost::noncopyable >
    ("PropertyBase", no_init)
    .def("checkFitsTo", &PropertyBase::checkFitsTo)
    ;

  class_< Property< Real3D >, boost::noncopyable, bases< PropertyBase > >
    ("Real3DProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< Real3D >::getItem)
    .def("__setitem__", &Property< Real3D >::setItem)
    ;
  
  class_< Property< int >, boost::noncopyable, bases< PropertyBase > >
    ("IntegerProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< int >::getItem)
    .def("__setitem__", &Property< int >::setItem)
    ;
  
  class_< Property< real >, boost::noncopyable, bases< PropertyBase > >
    ("RealProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Property< real >::getItem)
    .def("__setitem__", &Property< real >::setItem)
    ;
  
  class_< ArrayProperty< int >, boost::noncopyable, bases< PropertyBase > >
    ("IntegerArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &ArrayProperty< int >::getItem)
    .def("__setitem__", &ArrayProperty< int >::setItem)
    ;
  
  class_< ArrayProperty< real >, boost::noncopyable, bases< PropertyBase > >
    ("RealArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &ArrayProperty< real >::getItem)
    .def("__setitem__", &ArrayProperty< real >::setItem)
    ;
}

