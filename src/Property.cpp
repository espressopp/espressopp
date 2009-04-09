#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include "Property.hpp"
#include "Particle.hpp"

using namespace espresso::particles;

void espresso::registerPythonProperties()
{
  using namespace boost::python;

  class_< Property<Real3D> >("Real3DProperty",
                             init< boost::shared_ptr<Storage> >())
    .def("__getitem__", &Property<Real3D>::getItem)
    .def("__setitem__", &Property<Real3D>::setItem);

  class_< Property<int> >("IntegerProperty",
                          init< boost::shared_ptr<Storage> >())
    .def("__getitem__", &Property<int>::getItem)
    .def("__setitem__", &Property<int>::setItem);

  class_< Property<real> >("RealProperty",
                           init< boost::shared_ptr<Storage> >())
    .def("__getitem__", &Property<real>::getItem)
    .def("__setitem__", &Property<real>::setItem);

  class_< ArrayProperty<int> >("IntegerArrayProperty",
                               init< boost::shared_ptr<Storage>, size_t >())
    .def("__getitem__", &ArrayProperty<int>::getItem)
    .def("__setitem__", &ArrayProperty<int>::setItem);
  
  class_< ArrayProperty<real> >("RealArrayProperty",
                                init< boost::shared_ptr<Storage>, size_t >())
    .def("__getitem__", &ArrayProperty<real>::getItem)
    .def("__setitem__", &ArrayProperty<real>::setItem);
}

