#include <boost/python/implicit.hpp>
#include "python.hpp"
#include "RealND.hpp"

namespace espresso {
  struct real3D_pickle_suite : boost::python::pickle_suite {
    static
    python::tuple
    getinitargs(const RealND v){
      python::list ret;
      for(int i=0;  i<v.getDimension(); i++)
        ret.append( v[i] );
      return python::make_tuple( ret );
    }
  };
  
  void
  RealND::registerPython() {
    using namespace python;
    using namespace boost::python;
    class_< RealND >("RealND", init<>())
      .def(init< int >())
      .def(init< int, real >())
      .add_property("dimension", &RealND::getDimension, 
                                 &RealND::setDimension)
      .def("__getitem__", &RealND::getItem)
      .def("__setitem__", &RealND::setItem)
      .def("sqr", &RealND::sqr)
      .def("abs", &RealND::abs)
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
      //.def("cross", &RealND::cross)
      .def_pickle(real3D_pickle_suite())
      ;

  }
}
