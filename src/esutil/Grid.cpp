#include "Grid.hpp"
#include "python.hpp"

namespace espresso {
  namespace esutil {
    
    //BOOST_PYTHON_FUNCTION_OVERLOADS( testMethod1, Grid::testMethod, 1, 1)
    //BOOST_PYTHON_FUNCTION_OVERLOADS( testMethod2, Grid::testMethod, 1, 1)

    //void (Grid::*mapIndexToPositionInts)(int&, int&, int&, longint) const = &Grid::mapIndexToPosition;
    void (Grid::*mapIndexToPositionInt3D)(Int3D&, longint) const = &Grid::mapIndexToPosition;
    
    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void Grid::registerPython() {

      using namespace espresso::python;

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

