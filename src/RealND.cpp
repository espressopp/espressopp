/*
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
  
  This file is part of ESPResSo++.
  
  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/

#include <boost/python/implicit.hpp>
#include "python.hpp"
#include "RealND.hpp"

namespace espressopp {
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
