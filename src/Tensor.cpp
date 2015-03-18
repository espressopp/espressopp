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

#include <python.hpp>
#include <boost/python/implicit.hpp>
#include "Tensor.hpp"

namespace espressopp {
  struct tensor_pickle_suite : boost::python::pickle_suite {
    static
    python::tuple
    getinitargs(const Tensor v)
    { return python::make_tuple(v[0], v[1], v[2], v[3], v[4], v[5]); }
  };
  
  void
  Tensor::
  registerPython() {
    using namespace python;
    using namespace boost::python;
    class_< Tensor >("Tensor", init<>())
      .def(init< double, double, double, double, double, double >())
      .def("__getitem__", &Tensor::getItem)
      .def("__setitem__", &Tensor::setItem)
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
      .def_pickle(tensor_pickle_suite())
      ;

  }
}
