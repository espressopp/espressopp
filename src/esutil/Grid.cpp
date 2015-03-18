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

#include "python.hpp"
#include "Grid.hpp"

namespace espressopp {
  namespace esutil {
    
    //BOOST_PYTHON_FUNCTION_OVERLOADS( testMethod1, Grid::testMethod, 1, 1)
    //BOOST_PYTHON_FUNCTION_OVERLOADS( testMethod2, Grid::testMethod, 1, 1)

    //void (Grid::*mapIndexToPositionInts)(int&, int&, int&, longint) const = &Grid::mapIndexToPosition;
    void (Grid::*mapIndexToPositionInt3D)(Int3D&, longint) const = &Grid::mapIndexToPosition;
    
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void Grid::registerPython() {

      using namespace espressopp::python;

      class_< Grid >("esutil_Grid", init<>())
        .def(init< int, int, int>())
        .def(init<const Int3D&>())
        //.def("mapIndexToPosition", mapIndexToPositionInts)
        .def("mapIndexToPosition", mapIndexToPositionInt3D)
      /*
        .def("testMethod", 
              static_cast< void(Grid::*)(double) const>(&Grid::testMethod),
              testMethod1())

        .def("testMethod", 
              static_cast< void(Grid::*)(int) const>(&Grid::testMethod),
              testMethod2())
      */
      ;
    }
  }  
}

