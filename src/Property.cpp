#include "python.hpp"
#include "Property.hpp"
#include "Particle.hpp"

void espresso::registerPythonProperties()
{
  using namespace espresso::particles;
  using namespace espresso::python;

  class_< Real3DProperty, Real3DProperty::SelfPtr >
    ("Real3DProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &Real3DProperty::getItem)
    .def("__setitem__", &Real3DProperty::setItem);
  
  class_< IntegerProperty, IntegerProperty::SelfPtr >
    ("IntegerProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &IntegerProperty::getItem)
    .def("__setitem__", &IntegerProperty::setItem);
  
  class_< RealProperty, RealProperty::SelfPtr >
    ("RealProperty", init< Storage::SelfPtr >())
    .def("__getitem__", &RealProperty::getItem)
    .def("__setitem__", &RealProperty::setItem);
  
  class_< IntegerArrayProperty, IntegerArrayProperty::SelfPtr >
    ("IntegerArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &IntegerArrayProperty::getItem)
    .def("__setitem__", &IntegerArrayProperty::setItem);
  
  class_< RealArrayProperty, RealArrayProperty::SelfPtr >
    ("RealArrayProperty", init< Storage::SelfPtr, size_t >())
    .def("__getitem__", &RealArrayProperty::getItem)
    .def("__setitem__", &RealArrayProperty::setItem);
}

