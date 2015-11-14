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
#include "Quaternion.hpp"
#include "Real3D.hpp"

namespace espressopp {

  struct quaternion_pickle_suite : boost::python::pickle_suite {
    static
    python::tuple
    getinitargs(const Quaternion v)
    { return python::make_tuple(v.getReal(), v.getImag()); }
  };

  void
  Quaternion::
  registerPython() {
    using namespace python;
    using namespace boost::python;
    class_< Quaternion >("Quaternion", init<>())
      .def(init< double, double, double, double >())
      .def(init< double, Real3D >())
      .def(init< double >())
      .def(init< Real3D >())
      .add_property("real_part", &Quaternion::getReal, 
      	                 &Quaternion::setReal)
      .add_property("unreal_part", &Quaternion::getImag,
 	                   &Quaternion::setImag)
      .def("getReal", &Quaternion::getReal)
      .def("setReal", &Quaternion::setReal)
      .def("getImag", &Quaternion::getImag)
      .def("setImag", &Quaternion::setImag)
      .def("getImagItem", &Quaternion::getImagItem)
      .def("setImagItem", &Quaternion::setImagItem)
      .def("__setitem__", &Quaternion::setItem)
      .def("__getitem__", &Quaternion::getItem)
      .def("sqr", &Quaternion::sqr)
      .def("abs", &Quaternion::abs)
      .def("normalize", &Quaternion::normalize)
      .def("transpose", &Quaternion::transpose)
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
      .def_pickle(quaternion_pickle_suite())
      ;

  }
}
