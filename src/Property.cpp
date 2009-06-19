#include <boost/python.hpp>
#include "Property.hpp"
#include "Particle.hpp"

using namespace espresso::particles;

void espresso::registerPythonProperties()
{
  using namespace boost::python;

  class_< Real3DProperty >("Real3DProperty",
			   init< PStorage >())
    .def("__getitem__", &Real3DProperty::getItem)
    .def("__setitem__", &Real3DProperty::setItem);
  
  class_< IntegerProperty >("IntegerProperty",
			    init< PStorage >())
    .def("__getitem__", &IntegerProperty::getItem)
    .def("__setitem__", &IntegerProperty::setItem);
  
  class_< RealProperty >("RealProperty",
			 init< PStorage >())
    .def("__getitem__", &RealProperty::getItem)
    .def("__setitem__", &RealProperty::setItem);
  
  class_< ArrayProperty<int> >("IntegerArrayProperty",
                               init< PStorage, size_t >())
    .def("__getitem__", &IntegerArrayProperty::getItem)
    .def("__setitem__", &IntegerArrayProperty::setItem);
  
  class_< ArrayProperty<real> >("RealArrayProperty",
                                init< PStorage, size_t >())
    .def("__getitem__", &RealArrayProperty::getItem)
    .def("__setitem__", &RealArrayProperty::setItem);
}

